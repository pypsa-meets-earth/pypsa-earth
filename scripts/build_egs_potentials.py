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
        x_stepsize = transform[0]  # Width of a pixel in map units
        y_stepsize = abs(
            transform[4]
        )  # Height of a pixel in map units (absolute value since it's negative)
        logger.info(
            f"Longitude stepsize: {x_stepsize}, Latitude stepsize: {y_stepsize}"
        )

        assert x_stepsize == y_stepsize

    lats = np.arange(25, 50, y_stepsize) # US extent
    lons = np.arange(-125, -65, x_stepsize)

    area = np.mean(cell_areas(lats, lons, x_stepsize, y_stepsize))

    if mode == "egs":
        capacity_per_datapoint = area * 0.45 # assuming a 0.45MW/km2 footprint
    elif mode == "hs":
        capacity_per_datapoint = area * 0.024 # assuming a 0.024MW/km2 footprint
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

    gdf["heat_share"] = gdf["sales_heat"].div(gdf["sales_power"])

    # gdf["p_nom_max[MWe]"] = gdf["sales_power"] * plants_per_datapoint
    gdf["p_nom_max[MWe]"] = capacity_per_datapoint

    nodal_egs_potentials = pd.DataFrame(
        np.nan,
        index=regions.index,
        columns=["capex[USD/MWe]", "opex[USD/MWhe]", "p_nom_max[MWe]"],
    )

    config = snakemake.params["enhanced_geothermal"]

    regional_potentials = []

    for i, (name, geom) in tqdm(
        enumerate(regions.geometry.items()),
        desc="Matching EGS potentials to network regions",
        ascii=True,
    ):
        # For polygon geometries, we need to check which regions overlap with the polygon
        ss = gdf[
            gdf.geometry.intersects(geom)
        ]  # could be improved for partial overlaps
        if ss.empty:
            continue

        # Sort by capex first
        ss = (
            ss[["capex[USD/MWe]", "opex[USD/MWhe]", "p_nom_max[MWe]"]]
            .reset_index(drop=True)
            .sort_values(by="capex[USD/MWe]")
        )

        ss["agg_available_capacity[MWe]"] = ss["p_nom_max[MWe]"].cumsum()

        if ss.empty:
            continue

        # Instead of equal capex ranges, create bins with equal capacity
        if len(ss) > 1:
            total_capacity = ss["p_nom_max[MWe]"].sum()
            capacity_per_bin = total_capacity / config["supply_curve_steps"]

            # Create bins based on equal capacity distribution
            bin_edges = []
            current_capacity = 0
            target_capacity = capacity_per_bin

            # Create a dictionary to store data for each bin
            bin_data = []
            current_bin = {"capacity": 0, "capex_weighted_sum": 0}

            for _, row in ss.iterrows():
                current_bin["capacity"] += row["p_nom_max[MWe]"]
                current_bin["capex_weighted_sum"] += (
                    row["capex[USD/MWe]"] * row["p_nom_max[MWe]"]
                )

                # If we've reached or exceeded the target capacity for this bin
                if current_bin["capacity"] >= target_capacity:
                    # Calculate weighted average capex for the bin
                    weighted_capex = (
                        current_bin["capex_weighted_sum"] / current_bin["capacity"]
                    )
                    bin_data.append(
                        {
                            "level": myround(weighted_capex),
                            "p_nom_max[MWe]": current_bin["capacity"],
                            "opex[USD/MWhe]": row[
                                "opex[USD/MWhe]"
                            ],  # For simplicity, using the last opex
                        }
                    )

                    # Start a new bin
                    current_bin = {"capacity": 0, "capex_weighted_sum": 0}
                    target_capacity += capacity_per_bin

            # Add any remaining capacity to the last bin
            if current_bin["capacity"] > 0:
                weighted_capex = (
                    current_bin["capex_weighted_sum"] / current_bin["capacity"]
                )
                bin_data.append(
                    {
                        "level": myround(weighted_capex),
                        "p_nom_max[MWe]": current_bin["capacity"],
                        "opex[USD/MWhe]": ss["opex[USD/MWhe]"].iloc[
                            -1
                        ],  # Using the last opex
                    }
                )

            # Convert to DataFrame
            ss = pd.DataFrame(bin_data)
        else:
            ss["level"] = ss["capex[USD/MWe]"].iloc[0]

        ss = ss.dropna()

        # Convert to the desired format
        if len(ss) > 1:
            ss = ss.set_index("level")
        else:
            # Handle the single row case
            ss = ss.groupby("level")[["p_nom_max[MWe]", "opex[USD/MWhe]"]].agg(
                {"p_nom_max[MWe]": "sum", "opex[USD/MWhe]": "mean"}
            )

        ss.index = pd.MultiIndex.from_product(
            [[name], ss.index], names=["network_region", "capex[USD/MWe]"]
        )

        ss.loc[:, "p_nom_max[MWe]"] = ss.loc[:, "p_nom_max[MWe]"].iloc[0]

        regional_potentials.append(ss)

    regional_potentials = pd.concat(regional_potentials).dropna()

    regional_potentials = regional_potentials.reset_index(level="capex[USD/MWe]")

    regional_potentials["supply_curve_step"] = (
        regional_potentials.groupby("network_region").cumcount() + 1
    )

    regional_potentials = regional_potentials.set_index(
        "supply_curve_step", append=True
    )

    regional_potentials = regional_potentials.rename(
        columns={"capex[USD/MWe]": "capital_cost[USD/MWe]"}
    )

    regional_potentials.to_csv(snakemake.output.egs_potentials)
