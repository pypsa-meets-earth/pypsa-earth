# -*- coding: utf-8 -*-
"""Build EGS potentials for different depths and temperatures. Upper limits available
heat content by overlaying subsurface potential with geospatial demand data."""


import logging

logger = logging.getLogger(__name__)

import os

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from _helpers import configure_logging
from rasterio.transform import xy
from shapely.geometry import Point
from tqdm import tqdm


def get_raster_file(tif_file):

    with rasterio.open(tif_file) as src:

        band1 = src.read(1)
        transform = src.transform
        rows, cols = band1.shape

        row_indices, col_indices = np.meshgrid(
            np.arange(rows), np.arange(cols), indexing="ij"
        )
        row_indices = row_indices.flatten()
        col_indices = col_indices.flatten()
        band1 = band1.flatten()

        nodata = src.nodata
        valid_mask = band1 != nodata
        band1 = band1[valid_mask]
        row_indices = row_indices[valid_mask]
        col_indices = col_indices[valid_mask]

        xs, ys = xy(transform, row_indices, col_indices)

        df = pd.DataFrame({"value": band1, "x": xs, "y": ys})

        geometry = [Point(xy) for xy in zip(df["x"], df["y"])]

        return gpd.GeoDataFrame(df, geometry=geometry, crs=src.crs)


def myround(x):
    if x % 10 == 0:
        return int(x)
    elif x % 10 < 5:
        return int(x - x % 10)
    else:
        return int(x + 10 - x % 10)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_egs_potentials",
            simpl="",
            clusters=50,
        )

    configure_logging(snakemake)

    regions = gpd.read_file(snakemake.input.shapes).set_index("name")

    capex_gdf = (
        get_raster_file(snakemake.input.egs_capex)
        .set_index(["x", "y"])
        .rename(columns={"value": "capex"})
    )
    gen_gdf = (
        get_raster_file(snakemake.input.egs_gen)
        .set_index(["x", "y"])
        .rename(columns={"value": "gen"})
    )
    opex_gdf = (
        get_raster_file(snakemake.input.egs_opex)
        .set_index(["x", "y"])
        .rename(columns={"value": "opex"})
    )

    gdf = capex_gdf[["capex", "geometry"]].join(gen_gdf[["gen"]], how="inner")
    gdf = gdf.join(opex_gdf[["opex"]], how="inner")

    gdf["plant_capacity"] = gdf["gen"].div(8760)
    gdf["capex($/kWel)"] = gdf["capex"].div(gdf["plant_capacity"]).mul(1e3)
    gdf["plant_capacity(MWel)"] = gdf["plant_capacity"]
    gdf["opex[$/kWh]"] = gdf["opex"].mul(1e-2)

    gdf.rename(
        columns={
            "capex($/kWel)": "capex[$/kW]",
            "plant_capacity(MWel)": "available_capacity[MW]",
        },
        inplace=True,
    )

    gdf = gdf[["capex[$/kW]", "opex[$/kWh]", "available_capacity[MW]", "geometry"]]

    gdf = gdf.replace([np.inf, -np.inf], np.nan)

    nodal_egs_potentials = pd.DataFrame(
        np.nan,
        index=regions.index,
        columns=["capex[$/kW]", "opex[$/kWh]", "p_nom_max[MW]"],
    )

    config = snakemake.params["enhanced_geothermal"]

    regional_potentials = []

    for name, geom in tqdm(
        regions.geometry.items(),
        desc="Matching EGS potentials to network regions",
        ascii=True,
    ):

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
                # ss["capex[$/kW]"].min(),
                0,
                ss["capex[$/kW]"].max(),
                config["max_levels"] + 1,
            )
        )

        labels = bins.rolling(2).mean().dropna().tolist()

        if ss.empty:
            continue

        if len(ss) > 1:
            ss["level"] = pd.cut(
                ss["capex[$/kW]"], bins=bins, labels=labels, duplicates="drop"
            )
        else:
            ss["level"] = ss["capex[$/kW]"].iloc[0]

        ss = ss.dropna()

        ss = ss.groupby("level", observed=False)[
            ["available_capacity[MW]", "opex[$/kWh]"]
        ].agg({"available_capacity[MW]": "sum", "opex[$/kWh]": "mean"})
        ss.index = list(map(myround, ss.index))
        ss.index = pd.MultiIndex.from_product(
            [[name], ss.index], names=["network_region", "capex[$/kW]"]
        )

        regional_potentials.append(ss)

    regional_potentials = pd.concat(regional_potentials).dropna()
    dr = snakemake.params["costs"]["fill_values"]["discount rate"]

    def to_capital_cost(index, lt=25):
        return (index[0], index[1] * dr / (1 - (1 + dr) ** (-lt)))

    regional_potentials.index = regional_potentials.index.map(to_capital_cost)
    regional_potentials.index.names = ["network_region", "capital_cost[$/kW]"]

    regional_potentials.to_csv(snakemake.output.egs_potentials)
