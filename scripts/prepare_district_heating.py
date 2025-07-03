# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)

import geopandas as gpd
import numpy as np
import pandas as pd
from build_egs_potentials import tif_to_gdf
from build_industrial_heating_demand import process_regional_supply_curves

# from build_industrial_heating_demand import process_techno_economic_data
from tqdm import tqdm


def process_techno_economic_data(df):
    """
    Process the techno-economic data for different technologies.

    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing technology parameters

    Returns:
    --------
    pandas.DataFrame
        DataFrame with processed techno-economic parameters
    """
    # Define the standardized technology names mapping
    tech_mapping = {
        "directheat_100": "directheat100degC",
    }

    # Mapper for output types and their temperature bands
    tech_output_mapping = {
        "directheat100degC": {"heat": "50-80C"},
    }

    # Create a new DataFrame to store processed data
    result_data = {}
    result_data["geometry"] = df["geometry"]

    lifetimes = {"directheat": 30}

    # Process each technology type
    for tech_key, std_tech_name in tech_mapping.items():
        # Extract technology type and mode
        if "_egs" in tech_key:
            mode = "egs"
        elif "_hs" in tech_key:
            mode = "hs"
        else:
            mode = None

        # Calculate original LCOE from raw data
        temp = tech_key.split("_")[1]
        prefix = f"directheat_{temp}"
        lifetime = lifetimes["directheat"]

        # Get all capex, opex, and sales columns for this technology
        capex_cols = [
            col for col in df.columns if col.startswith(prefix) and "capex" in col
        ]
        opex_cols = [
            col for col in df.columns if col.startswith(prefix) and "opex" in col
        ]
        sales_cols = [
            col for col in df.columns if col.startswith(prefix) and "sales" in col
        ]

        if capex_cols and opex_cols and sales_cols:
            # Sum all capex, opex, and sales
            total_capex = sum(df[col] for col in capex_cols)
            total_opex = sum(df[col] for col in opex_cols)
            total_sales = sum(df[col] for col in sales_cols)

            # Calculate CRF (Capital Recovery Factor)
            discount_rate = 0.07
            crf = (
                discount_rate
                * (1 + discount_rate) ** lifetime
                / ((1 + discount_rate) ** lifetime - 1)
            )

            # Calculate original LCOE
            original_lcoe = (total_capex * crf + total_opex) / (total_sales * 8760)

            # Store the original LCOE
            # result_data[(std_tech_name, 'original_lcoe[USD/MWh]')] = original_lcoe * 1e6

        temp = tech_key.split("_")[1]

        # Column prefixes
        prefix = f"directheat_{temp}"

        # For directheat, sales, capex and opex are at the end of column names
        sales_col = f"{prefix}_sales"

        if sales_col not in df.columns:
            error_msg = f"Sales column {sales_col} not found for technology {tech_key}"
            print(error_msg)
            raise ValueError(error_msg)

        # Get the appropriate lifetime for this technology
        lifetime = lifetimes["directheat"]

        # Calculate total lifetime output
        lifetime_output = df[sales_col] * 8760 * lifetime
        result_data[(std_tech_name, "total_output[MWh]")] = lifetime_output

        # Get CAPEX
        capex_col = f"{prefix}_capex"
        if capex_col not in df.columns:
            error_msg = f"CAPEX column {capex_col} not found for technology {tech_key}"
            print(error_msg)
            raise ValueError(error_msg)
        result_data[(std_tech_name, "capex[USD/MW]")] = (
            df[capex_col].div(df[sales_col]).mul(1e6)
        )

        # Get OPEX
        opex_cols = [
            col for col in df.columns if col.startswith(prefix) and "opex" in col
        ]
        if not opex_cols:
            error_msg = (
                f"No OPEX columns found for technology {tech_key} with prefix {prefix}"
            )
            print(error_msg)
            raise ValueError(error_msg)
        # Sum all OPEX columns instead of just using the first one
        total_opex = sum(df[col] for col in opex_cols)
        result_data[(std_tech_name, "opex[USD/MWh]")] = total_opex.div(
            lifetime_output
        ).mul(1e6)

        # Only set the heat share where opex is not NaN
        heat_share_key = (
            std_tech_name,
            f"{tech_output_mapping[std_tech_name]['heat']}_share",
        )
        mask = ~total_opex.isna()
        result_data[heat_share_key] = pd.Series(np.nan, index=total_opex.index)
        result_data[heat_share_key].loc[mask] = 1.0

    # Create DataFrame with MultiIndex columns
    return pd.DataFrame(result_data)


if __name__ == "__main__":

    regions = gpd.read_file(snakemake.input.regions).set_index("name")
    district_gdf = gpd.read_file(snakemake.input.demand_data)

    # district_gdf["heating_demand_mwh"] = district_gdf["heating_by_pop"] * 2.93071e-7
    # district_gdf["avg_heat_mw"] = district_gdf["heating_demand_mwh"] / 8760

    district_gdf["avg_heat_mw"] = district_gdf["mwh"] / 8760

    # this relates to the piping-cost incured in each region
    district_gdf["boreholes_per_sqkm"] = district_gdf["avg_heat_mw"] / 10

    # we assume that a typical square kilometer of residential area that
    # can be supplied by district heating needs around 20km or piping
    # https://ojs.library.queensu.ca/index.php/cpp/article/download/13406/9382/31196
    # [ validated against Google Maps for DH-eligible area in Chicago ]

    # piping-cost [USD/m] * 20km / boreholes_per_sqkm = network-cost per borehole (i.e. network-cost per 10MW of borehole heat)

    # for piping-cost, we go with numbers from https://www.npro.energy/main/en/help/technology-costs
    # assuming Urban / Paved Surfaces
    # should be around 3250 Euro/m for a 0.30m diameter pipe
    # 1.15 is the conversion factor from Euro to USD
    piping_cost_per_m = 3250 * 1.15  # USD/m

    district_gdf["network_cost_per_mw"] = (
        piping_cost_per_m * 20_000 / district_gdf["boreholes_per_sqkm"] / 10
    )

    tif_files = {
        name: fn for name, fn in snakemake.input.items() if fn.endswith(".tif")
    }
    file_name_transformer = lambda x: "-".join(str(x).split("/")[-2:]).replace(
        ".tif", ""
    )

    gdf = tif_to_gdf(tif_files.values(), name_transformer=file_name_transformer)
    gdf = gdf.rename(
        columns={file_name_transformer(item): key for key, item in tif_files.items()}
    )

    gdf = process_techno_economic_data(gdf)

    gdf.rename(columns={"geometry": ("geometry", "")}, inplace=True)
    gdf.columns = pd.MultiIndex.from_tuples(gdf.columns)

    tech_levels = gdf.columns.get_level_values(0).unique()

    new_parts = [gdf[["geometry"]]]

    idx = pd.IndexSlice

    for tech in tech_levels:
        if tech == "geometry":
            continue

        new_parts.append(gdf.loc[:, idx[tech, :]].dropna())

    gdf = pd.concat(new_parts, axis=1)

    for tech in gdf.columns.get_level_values(0).unique()[::-1]:
        # Skip geometry column
        if tech == "geometry":
            continue

        # Get all columns for this technology
        idx = pd.IndexSlice
        tech_data = gdf.loc[:, idx[tech, :, :]].dropna()

        if tech_data.empty:
            continue

        # Calculate statistics, ignoring NaN values
        tech_averages = tech_data.mean(skipna=True)
        tech_medians = tech_data.median(skipna=True)
        tech_mins = tech_data.min(skipna=True)
        tech_maxs = tech_data.max(skipna=True)

        # Define discount rate and lifetimes for LCOE calculation
        discount_rate = 0.07
        lifetime = 30.0
        lcoe_values = {}

        tech_type = "directheat"

        # Find capex and opex columns for this technology
        capex_cols = [
            col
            for col in tech_data.columns.get_level_values(1)
            if "capex[USD/MW" in col
        ]
        opex_cols = [
            col
            for col in tech_data.columns.get_level_values(1)
            if "opex[USD/MWh" in col
        ]
        shares_cols = [
            col for col in tech_data.columns.get_level_values(1) if "share" in col
        ]

        # Get the first capex and opex columns (they should be for the same output type)
        assert len(capex_cols) == 1
        assert len(opex_cols) == 1

        capex_col = capex_cols[0]
        opex_col = opex_cols[0]

        # Extract the output type from the column name
        if not shares_cols:
            output_type = [tech]
        else:
            output_type = [col.split("_")[0] for col in shares_cols]

        # Calculate CRF (Capital Recovery Factor)
        crf = (
            discount_rate
            * (1 + discount_rate) ** lifetime
            / ((1 + discount_rate) ** lifetime - 1)
        )

        # Calculate LCOE for each row: LCOE = (CAPEX * CRF + OPEX)
        if shares_cols:
            shares = tech_data.loc[:, idx[:, shares_cols]].sum(axis=1)
            assert np.allclose(shares, 1)
        else:
            shares = 1
        lcoe = (
            tech_data.loc[:, idx[:, capex_col]].values.flatten() * crf / 8760
            + tech_data.loc[:, idx[:, opex_col]].values.flatten()
        )

        lcoe = pd.Series(lcoe, index=tech_data.index)

        # Add LCOE values to the gdf dataframe
        gdf.loc[lcoe.index, idx[tech, f"lcoe[USD/MWh]"]] = lcoe

    regional_supplies = list()

    # Convert string geometries to shapely objects
    if isinstance(gdf["geometry"].iloc[0], str):
        from shapely import wkt

        gdf["geometry"] = gdf["geometry"].apply(lambda x: wkt.loads(x))
    gdf = gpd.GeoDataFrame(gdf, geometry="geometry", crs="EPSG:4326")

    regional_supply_shapes = pd.Series(index=regions.index)

    total_results = []
    total_clusters = []

    for region, geometry in tqdm(regions["geometry"].items()):

        regional_supply = []

        ss = district_gdf.loc[district_gdf["geometry"].within(geometry)]

        if ss.empty:
            continue

        techs = [
            # "pwr_residheat80degC_egs",
            # "pwr_residheat80degC_hs",
            "directheat100degC",
        ]

        for index, row in ss[
            ["avg_heat_mw", "network_cost_per_mw", "geometry"]
        ].iterrows():

            query_point = row["geometry"]

            buffer_distance = 0.1  # in degrees
            buffered_point = query_point.buffer(buffer_distance)

            mask = gdf.geometry.intersects(buffered_point)
            if not sum(mask):
                continue

            geothermal_subset = gdf.loc[mask].iloc[0]

            if len(geothermal_subset) == 0 or geothermal_subset.isna().any():
                # no geothermal supply in this region
                continue

            geothermal_subset.index = geothermal_subset.index.get_level_values(1)

            supply_curve_step = pd.Series(
                {
                    "heat_demand[MW]": row["avg_heat_mw"],
                    "capex[USD/MW]": geothermal_subset.loc["capex[USD/MW]"]
                    + row["network_cost_per_mw"],
                    "opex[USD/MWh]": geothermal_subset.loc["opex[USD/MWh]"],
                }
            )
            regional_supply.append(supply_curve_step)

        if len(regional_supply) == 0:
            continue

        regional_supply = pd.concat(regional_supply, axis=1).T
        regional_supply["region"] = region
        regional_supply["tech"] = "directheat100degC"
        regional_supply[["lcoe[USD/MWh]", "total_output[MWh]"]] = np.nan

        regional_supply = process_regional_supply_curves(regional_supply)
        total_results.append(regional_supply)

    try:
        total_results = pd.concat(total_results)[
            ["heat_demand[MW]", "capex[USD/MW]", "opex[USD/MWh]"]
        ]
    except ValueError:
        raise ValueError("Total results is empty")
        total_results = pd.DataFrame(
            columns=["heat_demand[MW]", "capex[USD/MW]", "opex[USD/MWh]"]
        )

    total_results.to_csv(snakemake.output["district_heating_geothermal_supply_curves"])
