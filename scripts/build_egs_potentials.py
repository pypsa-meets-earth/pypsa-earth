# -*- coding: utf-8 -*-
"""Build EGS potentials for different depths and temperatures. Upper limits available
heat content by overlaying subsurface potential with geospatial demand data."""


import logging

logger = logging.getLogger(__name__)

import os
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from _helpers import configure_logging
from rasterio.transform import xy
from shapely.geometry import Point, box
from tqdm import tqdm

name_transformer = lambda x: Path(str(x).split("/")[-1]).name.replace(".tif", "")


def tif_to_gdf(tif_files, name_transformer=name_transformer):
    """
    Convert a list of .tif files into a GeoDataFrame where each pixel is represented as a square polygon.

    Parameters:
    -----------
    tif_files : list
        List of paths to .tif files

    Returns:
    --------
    geopandas.GeoDataFrame
        GeoDataFrame with pixel values and geometries
    """

    # Dictionary to store values for each geometry
    geometry_data = {}

    for tif_file in tqdm(tif_files, desc="Joining tif files to GeoDataFrame"):
        with rasterio.open(tif_file) as src:
            data = src.read(1)  # Read the first band
            transform = src.transform

            # Get the file name without path and without .tif extension
            source_name = name_transformer(tif_file)

            # Iterate through each pixel
            for row in range(data.shape[0]):
                for col in range(data.shape[1]):
                    value = data[row, col]

                    # Skip nodata values
                    if value == src.nodata or (src.nodata is None and value == 0):
                        continue

                    # Get the pixel's coordinates
                    x1, y1 = transform * (col, row)
                    x2, y2 = transform * (col + 1, row + 1)

                    # Create a box geometry for the pixel
                    geometry = box(min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2))

                    # Convert geometry to WKT for dictionary key
                    geom_wkt = geometry.wkt

                    # Initialize entry if this geometry hasn't been seen before
                    if geom_wkt not in geometry_data:
                        geometry_data[geom_wkt] = {"geometry": geometry}

                    # Add the value for this source
                    geometry_data[geom_wkt][source_name] = value

    # Convert dictionary to DataFrame
    df = pd.DataFrame.from_dict(geometry_data, orient="index")

    # Extract geometry column
    geometries = df.pop("geometry")

    # Create the GeoDataFrame
    gdf = gpd.GeoDataFrame(df, geometry=geometries, crs=src.crs)
    gdf.index = range(len(gdf))

    return gdf


def myround(x):
    if x % 10 == 0:
        return int(x)
    elif x % 10 < 5:
        return int(x - x % 10)
    else:
        return int(x + 10 - x % 10)


def cell_areas(lat, lon, x_stepsize, y_stepsize):
    """
    Approximate grid cell areas for lon/lat grid in sqkm.
    
    lat, lon: arrays of cell center coordinates (in degrees)
    x_stepsize, y_stepsize: grid resolution (in degrees)
    """
    R = 6371.0  # Earth radius in km

    # convert to radians
    lat_rad = np.radians(lat)
    dlon = np.radians(x_stepsize)
    dlat = np.radians(y_stepsize)

    # area formula on sphere: A = R^2 * dlon * (sin(lat+Δ/2) - sin(lat-Δ/2))
    area = (R**2) * dlon * (
        np.sin(lat_rad + dlat/2) - np.sin(lat_rad - dlat/2)
    )
    return area  # in sqkm


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_egs_potentials",
            simpl="",
            clusters=50,
        )

    configure_logging(snakemake)
    mode = snakemake.wildcards.mode

    regions = gpd.read_file(snakemake.input.shapes).set_index("name")

    tif_files = {
        name: fn for name, fn in snakemake.input.items() if fn.endswith(".tif")
    }

    gdf = tif_to_gdf(tif_files.values())

    with rasterio.open(list(tif_files.values())[0]) as src:
        transform = src.transform
        x_stepsize = transform[0]
        y_stepsize = abs(
            transform[4]
        )

        assert x_stepsize == y_stepsize

    if mode == "egs":
        capacity_density = 0.45 # assuming a 0.45MW/km2 footprint
    elif mode == "hs":
        capacity_density = 0.024 # assuming a 0.024MW/km2 footprint
    else:
        raise ValueError(f"Invalid mode: {mode}")

    gdf = gdf.rename(
        columns={
            item.split("/")[-1].replace(".tif", ""): key
            for key, item in tif_files.items()
        }
    )

    # The following defines the power capacity is the measure of the
    # the plants capacity in the model, and defines the heat capacity
    # as the share of the power capacity that is used on top for heating.
    # (i.e this will be a multilink that dispatches more than one unit of
    # energy across links)
    gdf["capex[USD/MWe]"] = (
        (gdf["capex_power"] + gdf["capex_heat"] + gdf["capex_subsurf"])
        .div(gdf["sales_power"])
        .mul(1e6)
    )

    gdf["total_output[MWhe]"] = (
        gdf["sales_power"] * 8760 * 25
    )  # total energy output over lifetime
    gdf["opex[USD/MWhe]"] = (
        (gdf["opex_power"] + gdf["opex_heat"] + gdf["opex_subsurf"])
        .div(gdf["total_output[MWhe]"])
        .mul(1e6)
    )

    config = snakemake.params["enhanced_geothermal"]
    n_steps = config["supply_curve_steps"]

    regional_potentials = []

    for i, (name, geom) in tqdm(
        enumerate(regions.geometry.items()),
        desc="Matching EGS potentials to network regions",
        ascii=True,
    ):

        ss = gdf[
            gdf.geometry.intersects(geom)
        ]

        if ss.empty:
            continue

        ss = (
            ss[["capex[USD/MWe]", "opex[USD/MWhe]"]]
            .reset_index(drop=True)
            .sort_values(by="capex[USD/MWe]")
        )

        centroid = geom.centroid
        area = cell_areas(centroid.y, centroid.x, x_stepsize, y_stepsize)
        ss["p_nom_max[MWe]"] = capacity_density * area

        binned = ss.groupby(np.arange(len(ss)) // (len(ss) / n_steps)).agg({
            "capex[USD/MWe]": "mean",
            "opex[USD/MWhe]": "mean",
            "p_nom_max[MWe]": "sum",
        })
        binned.index = pd.MultiIndex.from_product(
            [[name], range(len(binned))], names=["network_region", "supply_curve_step"]
        )

        regional_potentials.append(binned)

    regional_potentials = pd.concat(regional_potentials)
    assert regional_potentials.notna().all().all(), "There are NaNs in the regional potentials"

    regional_potentials = regional_potentials.rename(
        columns={"capex[USD/MWe]": "capital_cost[USD/MWe]"}
    )

    total_capacity = regional_potentials['p_nom_max[MWe]'].sum() * 1e-6
    logger.info(f"Found {total_capacity:.3f} TWel {mode} potential. Saving to {snakemake.output.egs_potentials}.")
    regional_potentials.to_csv(snakemake.output.egs_potentials)
