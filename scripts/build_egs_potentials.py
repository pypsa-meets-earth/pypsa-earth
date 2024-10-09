# -*- coding: utf-8 -*-
"""
Build EGS potentials for different depths and temperatures.

Upper limits available heat content by overlaying subsurface potential
with geospatial demand data.
"""


import os

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm


def get_demands(df, x, y):
    raise NotImplementedError("implement me")


def get_demand_geometries(demand):
    raise NotImplementedError("implement me")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "build_egs_potentials",
            simpl="",
            clusters=50,
        )

    regions = gpd.read_file(snakemake.input.regions_onshore).set_index("name")
    egs_potentials = pd.read_csv(snakemake.input["egs_potential"], index_col=[0, 1, 2])

    gdf = gpd.GeoDataFrame(
        egs_potentials,
        geometry=gpd.points_from_xy(
            egs_potentials.index.get_level_values("lon"),
            egs_potentials.index.get_level_values("lat"),
        ),
    ).set_crs(epsg=4326)

    nodal_egs_potentials = pd.DataFrame(
        np.nan,
        index=regions.index,
        columns=["capex[$/kW]", "opex[$/kWh]", "p_nom_max[MW]"],
    )

    config = snakemake.params["enhanced_geothermal"]

    regional_potentials = []

    for name, geom in regions.geometry.items():

        ss = gdf.loc[gdf.geometry.within(geom)]
        if ss.empty:
            continue

        ss = (
            ss[["capex[$/kW]", "opex[$/kWh]", "available_capacity[MW]"]]
            .reset_index(drop=True)
            .sort_values(by="capex[$/kW]")
        )

        ss["agg_available_capacity[MW]"] = ss["available_capacity[MW]"].cumsum()

        bins = pd.Series(
            np.linspace(
                ss["capex[$/kW]"].min(),
                ss["capex[$/kW]"].max(),
                config["max_levels"] + 1,
            )
        )

        labels = bins.rolling(2).mean().dropna().tolist()

        ss["level"] = pd.cut(ss["capex[$/kW]"], bins=bins, labels=labels)
        ss = ss.dropna()

        ss = ss.groupby("level")[["available_capacity[MW]", "opex[$/kWh]"]].agg(
            {"available_capacity[MW]": "sum", "opex[$/kWh]": "mean"}
        )
        ss.index = pd.MultiIndex.from_product(
            [[name], ss.index], names=["network_region", "capex[$/kW]"]
        )

        regional_potentials.append(ss)

    pd.concat(regional_potentials).dropna().to_csv(snakemake.output.egs_potentials)
