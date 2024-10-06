# -*- coding: utf-8 -*-
"""Build EGS potentials for different depths and temperatures. Upper limits available
heat content by overlaying subsurface potential with geospatial demand data."""


import os

import numpy as np
import xarray as xr
import pandas as pd
from tqdm import tqdm
import geopandas as gpd



def get_demands(df, x, y):
    raise NotImplementedError('implement me')


def get_demand_geometries(demand):
    raise NotImplementedError('implement me')


if __name__ == "__main__":
    if "snakemake" not in globals():
        from helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "build_egs_potentials",
            simpl="",
            clusters=38,
        )
        sets_path_to_root("pypsa-earth-sec")


    regions = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    demands = ["AC", "industrial_hot", "industrial_medium"]
    # demand_array = xr.open_dataarray(snakemake.input.demand_array)

    # egs_potentials = pd.read_csv('egs_dummy_data.csv', index_col=[0,1,2])
    egs_potentials = pd.read_csv(snakemake.input["egs_potential"], index_col=[0,1,2])
    gdf = gpd.GeoDataFrame(
        egs_potentials,
        geometry=gpd.points_from_xy(
            egs_potentials.index.get_level_values('lon'),
            egs_potentials.index.get_level_values('lat')
            )).set_crs(epsg=4326)

    node_egs_capex = pd.DataFrame(index=regions.index, columns=demands)

    config = snakemake.params["enhanced_geothermal"]

    regional_potentials = []


    for name, geom in regions.geometry.items():

        ss = gdf.loc[gdf.geometry.within(geom)]
        if ss.empty:
            continue

        ss = (
            ss[['capex', 'opex']]
            .reset_index(drop=True)
            .sort_values(by='capex')
        )

        ss['potential'] = ss.index
        ss['agg_potential'] = ss['potential'].cumsum()

        bins = pd.Series(np.linspace(ss['capex'].min(), ss['capex'].max(), config['max_levels']+1))
        labels = bins.rolling(2).mean().dropna().tolist()

        ss['level'] = pd.cut(ss['capex'], bins=bins, labels=labels)

        ss = (
            ss
            .groupby('level')[['potential', 'opex']]
            .agg({'potential': 'sum', 'opex': 'mean'})
        )
        ss.index = pd.MultiIndex.from_product([[name], ss.index], names=['network_region', 'capex'])

        regional_potentials.append(ss) 

        # to be included once heat usage data is available
        """
        geom = regions.loc[region, 'geometry']
        demand_geoms = get_demand_geometries(demand_array, geom)

        rings = get_circles(
            geom,
            config['demand_steps'],
            config['demand_max_distance'],
        )

        for ring in rings:

            for demand in demands:
                demand = get_demands(demand_array, geom, demand_geoms)
        """

    pd.concat(regional_potentials).dropna().to_csv(snakemake.output.egs_potentials)
