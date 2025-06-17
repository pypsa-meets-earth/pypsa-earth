# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd
import geopandas as gpd

from tqdm import tqdm

from build_egs_potentials import tif_to_gdf
from build_industrial_heating_demand import process_techno_economic_data

if __name__ == "__main__":

    regions = gpd.read_file(
        snakemake.input.regions
    )
    district_gdf = gpd.read_file(
        snakemake.input.demand_data
    )
    district_gdf['heating_demand_mwh'] = district_gdf['heating_by_pop'] * 2.93071e-7
    district_gdf['avg_heat_mw'] = district_gdf['heating_demand_mwh'] / 8760
    
    tif_files = {
        name: fn for name, fn in snakemake.input.items() if fn.endswith(".tif")
    }
    file_name_transformer = lambda x: "-".join(str(x).split("/")[-2:]).replace(
        ".tif", ""
    )

    # gdf = tif_to_gdf(tif_files.values(), name_transformer=file_name_transformer)
    # gdf = gdf.rename(
    #     columns={file_name_transformer(item): key for key, item in tif_files.items()}
    # )

    # gdf = process_techno_economic_data(gdf)
    # print(gdf.head())
    # gdf.to_csv('hold.csv')

    gdf = pd.read_csv('hold.csv', index_col=0)

    regional_supplies = list()

    # Convert string geometries to shapely objects
    if isinstance(gdf["geometry"].iloc[0], str):
        from shapely import wkt

        gdf["geometry"] = gdf["geometry"].apply(lambda x: wkt.loads(x))
    gdf = gpd.GeoDataFrame(gdf, geometry="geometry", crs="EPSG:4326")
    print(gdf.head())

    regional_supply_shapes = pd.Series(index=regions.index)

    total_results = []
    total_clusters = []

    final_demands = pd.DataFrame(
        index=regions.index,
        columns=["district_heating_demand[MW]"],
    )

    for region, geometry in tqdm(regions["geometry"].items()):

        regional_supply = []

        ss = district_gdf.loc[district_gdf["geometry"].within(geometry)]

        print(ss.head())
    
        techs = [
            # "pwr_residheat80degC_egs",
            # "pwr_residheat80degC_hs",
            "directheat100degC",
        ]

        for index, row in ss[['avg_heat_mw', 'geometry']].iterrows():
        
            query_point = row['geometry']

            buffer_distance = 0.1  # in degrees
            buffered_point = query_point.buffer(buffer_distance)

            print(gdf.head())
            geothermal_subset = gdf.loc[gdf.geometry.intersects(buffered_point)]
            print(geothermal_subset.head())

            if len(geothermal_subset) == 0:
                continue
            else:
                geothermal_subset = geothermal_subset.iloc[0]

            final_demands.loc[region, "district_heating_demand[MW]"] += row['avg_heat_mw']

            idx = pd.IndexSlice

            district_supply = (
                geothermal_subset.loc[idx[techs, :]]
                .dropna()
                .unstack()
                .replace(np.nan, 0)
            )

            print(district_supply.head())
            if (
                district_supply.empty
            ):
                continue

            continue

            cluster_supply["band_lcoe"] = cluster_supply["lcoe[USD/MWh]"].mul(
                1 - cluster_supply[f"{band}_share"]
            )

            cluster_supply = cluster_supply.sort_values(
                by=["band_lcoe", "lcoe[USD/MWh]"]
            ).iloc[[0]]
            cluster_supply.drop(columns=["band_lcoe"], inplace=True)

            cluster_supply.loc[cluster_supply.index[0], ["heat_demand[MW]"]] = (
                cluster_size
            )
            cluster_supply.loc[cluster_supply.index[0], "capex[USD/MW]"] = (
                cluster_supply.loc[cluster_supply.index[0], "capex[USD/MW]"]
                * cluster_size
                + piping_cost * cluster_sp
            ) / cluster_size

            regional_supply.append(cluster_supply)

            tech_data = tech.loc[tech["geometry"].within(geometry)]
            print(tech_data.head())
            regional_supply.append(tech_data["avg_heat_mw"].sum())
        
        regional_supply_shapes.loc[region] = regional_supply

    import sys
    sys.exit()
