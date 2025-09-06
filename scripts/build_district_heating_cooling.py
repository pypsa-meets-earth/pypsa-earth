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

    district_gdf["avg_heat_mw"] = district_gdf["frac_htg"] / 8760
    district_gdf["avg_cooling_mw"] = district_gdf["frac_clg"] / 8760

    # remove regions with less than 2 MW of heat or 1.5 MW of cooling demand, as these are too small to be viable
    # these are small numbers because the model can decide against installing geothermal DH/DC itself
    district_gdf = district_gdf[(district_gdf["avg_heat_mw"] > 2) | (district_gdf["avg_cooling_mw"] > 1.5)]

    # References for heat-density classification and conversions used below:
    # - Heat Roadmap Europe / STRATEGO (WP2 Country Heat Roadmaps): DH suitability classes
    #   0–30, 30–100, 100–300, >300 TJ/km²·yr. See Connolly et al., *Energy* (2014) and
    #   Paardekooper et al., Heat Roadmap Europe 4 reports (2016–2018).
    # - Unit conversion: 1 MW·yr = 31.536 TJ  ⇒  MW/sqkm = (TJ/sqkm·yr) / 31.536.
    # - Linear heat-density viability rule-of-thumb /leq 1.5 MWh/(m·yr) for DH networks:
    #   Persson & Werner, *Applied Energy* 88 (2011); Danish Energy Agency,
    #   *Technology Data for District Heating* (latest edition).

    # --- Research-based heat-density tiers (HRE/STRATEGO -> MW/km²)
    # Original classes: 0–30, 30–100, 100–300, >300 TJ/km²·yr
    # Convert TJ/yr -> MW using 1 MW·yr = 31.536 TJ
    bounds_mw = [0.0, 30/31.536, 100/31.536, 300/31.536, np.inf]
    labels = ["low", "med", "high", "very_high"]

    district_gdf["density_tier"] = pd.cut(
    district_gdf["avg_heat_mw"],
    bins=bounds_mw, labels=labels, right=True, include_lowest=True
    )

    # for piping-cost, we go with numbers from https://www.npro.energy/main/en/help/technology-costs
    # and assume a mixture of pipe diameters (with branches subject to smaller pipe diameters)
    # --- Installed $/m: cheaper outside cores, higher in dense urban ROWs
    USD_PER_M = {
        "very_high": 2200.0,  # dense downtown cores
        "high":      1400.0,  # mixed mid-rise
        "med":       900.0,  # urban fringe/smaller towns
        "low":       600.0,  # single-family/greenfield-like installs
    }
    district_gdf["piping_cost_per_m"] = (
        district_gdf["density_tier"].map(USD_PER_M).astype(float)
    )

    TRENCH_M_PER_KM2 = {"low": 8_000, "med": 12_000, "high": 15_000, "very_high": 20_000}
    district_gdf["trench_m_per_km2"] = district_gdf["density_tier"].map(TRENCH_M_PER_KM2).astype(float)

    district_gdf["heating_network_cost_per_mw"] = (
        district_gdf["piping_cost_per_m"] * district_gdf["trench_m_per_km2"] / district_gdf["avg_heat_mw"]
    )

    district_gdf["cooling_network_cost_per_mw"] = (
        district_gdf["piping_cost_per_m"] * district_gdf["trench_m_per_km2"] / district_gdf["avg_cooling_mw"]
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

    # Convert string geometries to shapely objects
    if isinstance(gdf["geometry"].iloc[0], str):
        from shapely import wkt

        gdf["geometry"] = gdf["geometry"].apply(lambda x: wkt.loads(x))
    gdf = gpd.GeoDataFrame(gdf, geometry="geometry", crs="EPSG:4326")

    regional_supply_shapes = pd.Series(index=regions.index)

    heating_total_results = []
    cooling_total_results = []

    for region, geometry in tqdm(regions["geometry"].items()):

        heating_regional_supply = []
        cooling_regional_supply = []

        ss = district_gdf.loc[district_gdf["geometry"].within(geometry)]

        if ss.empty:
            continue

        techs = [
            # "pwr_residheat80degC_egs",
            # "pwr_residheat80degC_hs",
            "directheat100degC",
        ]

        for index, row in ss[
            [
                "avg_heat_mw",
                "heating_network_cost_per_mw",
                "avg_cooling_mw",
                "cooling_network_cost_per_mw",
                "geometry",
            ]
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

            heating_supply_curve_step = pd.Series(
                {
                    "heat_demand[MW]": row["avg_heat_mw"],
                    "capex[USD/MW]": geothermal_subset.loc["capex[USD/MW]"]
                    + row["heating_network_cost_per_mw"],
                    "opex[USD/MWh]": geothermal_subset.loc["opex[USD/MWh]"],
                }
            )
            heating_regional_supply.append(heating_supply_curve_step)

            cooling_supply_curve_step = pd.Series(
                {
                    "cooling_demand[MW]": row["avg_cooling_mw"],
                    "capex[USD/MW]": geothermal_subset.loc["capex[USD/MW]"]
                    + row["cooling_network_cost_per_mw"],
                    "opex[USD/MWh]": geothermal_subset.loc["opex[USD/MWh]"],
                }
            )
            cooling_regional_supply.append(cooling_supply_curve_step)

        if len(heating_regional_supply) == 0:
            continue

        heating_regional_supply = pd.concat(heating_regional_supply, axis=1).T
        heating_regional_supply["region"] = region
        heating_regional_supply["tech"] = "directheat100degC"
        heating_regional_supply[["total_output[MWh]"]] = np.nan

        heating_regional_supply = process_regional_supply_curves(
            heating_regional_supply,
            demand_column="heat_demand[MW]"
        )

        heating_total_results.append(heating_regional_supply)

        cooling_regional_supply = pd.concat(cooling_regional_supply, axis=1).T
        cooling_regional_supply["region"] = region
        cooling_regional_supply["tech"] = "directcooling100degC"
        cooling_regional_supply[["total_output[MWh]"]] = np.nan

        cooling_regional_supply = process_regional_supply_curves(
            cooling_regional_supply,
            demand_column="cooling_demand[MW]"
        )

        cooling_total_results.append(cooling_regional_supply)
    try:
        total_results_heating = pd.concat(heating_total_results)[
            ["heat_demand[MW]", "capex[USD/MW]", "opex[USD/MWh]"]
        ]
    except ValueError:
        total_results_heating = pd.DataFrame(
            columns=["heat_demand[MW]", "capex[USD/MW]", "opex[USD/MWh]"]
        )

    try:
        total_results_cooling = pd.concat(cooling_total_results)[
            ["cooling_demand[MW]", "capex[USD/MW]", "opex[USD/MWh]"]
        ]
    except ValueError:
        total_results_cooling = pd.DataFrame(
            columns=["cooling_demand[MW]", "capex[USD/MW]", "opex[USD/MWh]"]
        )

    total_results_heating.to_csv(snakemake.output["district_heating_geothermal_supply_curves"])

    # for district cooling, add cost of absorption chillers
    # data taken from https://www.energy.gov/sites/prod/files/2017/06/f35/CHP-Absorption%20Chiller-compliant.pdf
    cop = 0.72
    capex = 2_100_000.
    opex = 1.

    total_results_cooling["capex[USD/MW]"] = (total_results_cooling["capex[USD/MW]"] + capex) / cop
    total_results_cooling["opex[USD/MWh]"] = total_results_cooling["opex[USD/MWh]"] + opex

    # expressed as the MW that can be provided in terms of heat
    total_results_cooling["cooling_demand[MW]"] = total_results_cooling["cooling_demand[MW]"] / cop

    total_results_cooling.to_csv(snakemake.output["district_cooling_geothermal_supply_curves"])
