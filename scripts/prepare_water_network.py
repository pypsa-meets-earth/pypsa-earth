# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Prepare water network for Desalination and input to Electrolysis in prepare_sector_network.
"""

import logging

logger = logging.getLogger(__name__)

import io
import os
import zipfile
from pathlib import Path

import country_converter as coco
import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import requests
import xarray as xr
from _helpers import (
    BASE_DIR,
    content_retrieve,
    progress_retrieve,
    two_2_three_digits_country,
)
from geopandas import GeoSeries
from matplotlib.colors import LinearSegmentedColormap, to_rgba
from rasterio import sample
from rasterio.features import shapes
from rasterio.mask import mask
from rasterio.plot import show
from rasterio.warp import Resampling, calculate_default_transform, reproject
from shapely.geometry import LineString, shape
from shapely.ops import nearest_points


def download_gebco():
    """
    Downloads the GEBCO bathymetric data, extracts it to a local directory, and saves it to snakemake.output.
    """
    url = "https://www.bodc.ac.uk/data/open_download/gebco/gebco_2024_sub_ice_topo/zip/"
    gebco_file = "GEBCO_2024_sub_ice_topo.nc"

    zip_fn = Path(os.path.join(BASE_DIR, "gebco.zip"))
    to_fn = Path(os.path.join(BASE_DIR, "data/gebco"))

    # Path to the extracted file
    gebco_file_path = to_fn / gebco_file

    # Check if the target file already exists
    if gebco_file_path.exists():
        logger.info(
            f"Gebco data already exists at '{gebco_file_path}'. Skipping download."
        )
        return gebco_file_path

    logger.info(f"Downloading GEBCO data from '{url}'.")
    try:
        progress_retrieve(url, zip_fn)
        logger.info(f"Downloaded GEBCO zip file to '{zip_fn}'.")
    except Exception as e:
        logger.error(f"Failed to download GEBCO data: {e}")
        return None

    logger.info(f"Extracting GEBCO data to '{to_fn}'.")
    with zipfile.ZipFile(zip_fn) as z:
        z.extractall(to_fn)
        logger.info("Extraction complete.")

    zip_fn.unlink()  # Delete the temporary zip file
    logger.info(f"GEBCO data available in '{to_fn}'.")

    return gebco_file_path


def download_aqueduct():
    """
    Downloads the Aqueduct 4.0 Water Risk data, extracts it to a local directory,
    and returns it.
    """
    url = "https://files.wri.org/aqueduct/aqueduct-4-0-water-risk-data.zip"
    aqueduct_file = "Aqueduct40_waterrisk_download_Y2023M07D05/GDB/Aq40_Y2023D07M05.gdb"

    zip_fn = Path(os.path.join(BASE_DIR, "aqueduct.zip"))
    to_fn = Path(os.path.join(BASE_DIR, "data/aqueduct"))

    aqueduct_file_path = to_fn / aqueduct_file  # Full path to the extracted target file

    # Check if the target file already exists
    if aqueduct_file_path.exists():
        logger.info(
            f"Aqueduct data already exists at '{aqueduct_file_path}'. Skipping download."
        )
        return aqueduct_file_path

    logger.info(f"Downloading Aqueduct data from '{url}'.")
    try:
        progress_retrieve(url, zip_fn)
        logger.info(f"Downloaded Aqueduct zip file to '{zip_fn}'.")
    except Exception as e:
        logger.error(f"Failed to download Aqueduct data: {e}")
        return None

    logger.info(f"Extracting Aqueduct data to '{to_fn}'.")
    try:
        with zipfile.ZipFile(zip_fn) as z:
            z.extractall(to_fn)
            logger.info("Extraction complete.")
    except Exception as e:
        logger.error(f"Failed to extract Aqueduct data: {e}")
        return None

    zip_fn.unlink()  # Delete the temporary zip file
    logger.info(f"Aqueduct data available in '{to_fn}'.")

    return aqueduct_file_path


def download_shorelines():
    """
    Downloads the Global Self-consistent, Hierarchical, High-resolution Geography Database (GSHHG)
    shapefile archive, extracts the relevant shoreline shapefile, and saves it to snakemake.output.
    """
    url = (
        "https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-2.3.7.zip"
    )
    shoreline_file = "GSHHS_shp/i/GSHHS_i_L1.shp"  # Relative path inside the zip
    output_file = snakemake.output.shorelines

    # Save locations
    zip_fn = Path(os.path.join(BASE_DIR, "shorelines.zip"))
    to_fn = Path(os.path.join(BASE_DIR, "data/shorelines"))

    # Path to the extracted shapefile
    extracted_shoreline_path = to_fn / shoreline_file

    # Check if the target file already exists
    if extracted_shoreline_path.exists():
        logger.info(
            f"Shoreline data already exists at '{extracted_shoreline_path}'. Skipping download."
        )

        # Load the shapefile using GeoPandas
        shoreline_gdf = gpd.read_file(extracted_shoreline_path)

        # Save the extracted shapefile to snakemake.output
        shoreline_gdf.to_file(output_file)

        return shoreline_gdf

    logger.info(f"Downloading databundle from '{url}'.")
    progress_retrieve(url, zip_fn)

    logger.info(f"Extracting databundle.")
    # Extract the zip file into memory
    with zipfile.ZipFile(zip_fn) as z:
        # Check if the specific file exists in the archive
        if shoreline_file in z.namelist():
            logger.info(f"Extracting {shoreline_file}.")
            # Extract only the specific shapefile and its related files
            for file in z.namelist():
                if file.startswith(shoreline_file.rsplit("/", 1)[0]):
                    z.extract(file, to_fn)
            logger.info("Extraction complete.")
        else:
            logger.info(f"File {shoreline_file} not found in the archive.")
            return None

    zip_fn.unlink()
    logger.info(f"Shorelines data available in '{to_fn}'.")

    # Load the shapefile using GeoPandas
    shoreline_gdf = gpd.read_file(extracted_shoreline_path)

    # Save the extracted shapefile to snakemake.output
    shoreline_gdf.to_file(output_file)

    return shoreline_gdf


def clip_shorelines_country(shoreline_gdf, country_shapes):
    """
    Clips the shorelines to the boundaries of a given country with a buffer applied.

    Parameters:
    - shoreline_gdf (GeoDataFrame): GeoDataFrame containing shoreline geometries.
    - country_shapes (GeoDataFrame): GeoDataFrame containing country boundary geometries.

    Returns:
    - clipped_shoreline (GeoDataFrame): GeoDataFrame containing the clipped shoreline geometries.

    Functionality:
    - Ensures the country GeoDataFrame is in a projected CRS (e.g., UTM) for accurate buffering.
    - Applies a buffer of 3000 meters to the country boundaries.
    - Converts the country GeoDataFrame back to its original CRS (EPSG:4326).
    - Ensures CRS consistency between the shoreline and country GeoDataFrames.
    - Extracts the boundary (exterior) of the shoreline polygons.
    - Clips the shoreline boundary to the buffered country boundary using a spatial intersection.
    """
    country_gdf = country_shapes.copy()

    # Ensure the GeoDataFrame is in a projected CRS (e.g., UTM)
    if country_gdf.crs.is_geographic:
        country_gdf = country_gdf.to_crs(epsg=3857)  # Web Mercator (meters)

    # Apply buffer (3000 meters)
    country_gdf["geometry"] = country_gdf.buffer(3000)

    # Convert back to the original CRS if needed
    country_gdf = country_gdf.to_crs(epsg=4326)

    # Ensure CRS consistency
    if shoreline_gdf.crs != country_gdf.crs:
        shoreline_gdf = shoreline_gdf.to_crs(country_gdf.crs)

    # Extract the boundary (exterior) of the shoreline polygons
    shoreline_boundary = shoreline_gdf.boundary

    # Convert the boundary to a GeoDataFrame
    shoreline_boundary_gdf = gpd.GeoDataFrame(
        geometry=shoreline_boundary, crs=shoreline_gdf.crs
    )

    # Clip the shoreline boundary to the country boundary
    clipped_shoreline = gpd.overlay(
        shoreline_boundary_gdf, country_gdf, how="intersection"
    )

    return clipped_shoreline


def clip_shorelines_natura(natura_tiff_path, country_shapes, clipped_shoreline):
    """
    Clips the shorelines based on Natura 2000 raster data and a 4 km buffer areas.

    Parameters:
    - natura_tiff_path (str): Path to the Natura 2000 raster file.
    - country_shapes (GeoDataFrame): GeoDataFrame containing country boundary geometries.
    - clipped_shoreline (GeoDataFrame): GeoDataFrame containing previously clipped shoreline geometries.

    Returns:
    - clipped_shoreline_natura (GeoDataFrame): GeoDataFrame containing the shoreline geometries after removing buffered Natura 2000 areas.

    Functionality:
    - Ensures CRS consistency between the raster, country shapes, and shoreline GeoDataFrames.
    - Clips the raster using the country geometry and extracts positive value areas.
    - Converts positive raster areas into polygons and buffers them by 4 km.
    - Removes buffered Natura 2000 areas from the clipped shoreline using spatial difference.
    - Ensures the final GeoDataFrames are in EPSG:4326 for plotting.
    - Saves the resulting GeoDataFrames to files specified by `snakemake.output`.
    """
    # Ensure CRS of shoreline matches the raster and country_shapes
    with rasterio.open(natura_tiff_path) as src:
        if country_shapes.crs != src.crs:
            country_shapes = country_shapes.to_crs(src.crs)
        if clipped_shoreline.crs != src.crs:
            clipped_shoreline = clipped_shoreline.to_crs(src.crs)

        # Clip the raster using the country geometry
        country_geometry = [country_shapes.unary_union]
        clipped_raster, clipped_transform = mask(src, country_geometry, crop=True)

        # Extract positive values (greater than 0) as a mask
        positive_mask = clipped_raster[0] > 0  # Assuming single-band raster
        shapes_gen = shapes(positive_mask.astype(np.uint8), transform=clipped_transform)

        # Convert raster positive value areas into polygons
        positive_polygons = [shape(geom) for geom, val in shapes_gen if val == 1]

    # Create a GeoDataFrame from the positive polygons
    positive_area_gdf = gpd.GeoDataFrame({"geometry": positive_polygons}, crs=src.crs)

    # Switch to a projected CRS for buffering (e.g., EPSG:3857)
    projected_crs = "EPSG:3857"
    positive_area_gdf = positive_area_gdf.to_crs(projected_crs)
    clipped_shoreline = clipped_shoreline.to_crs(projected_crs)

    # Buffer the positive areas by 4 km (4000 meters)
    buffered_positive_area_gdf = positive_area_gdf.buffer(4000)

    # Drop empty or invalid geometries after buffering
    buffered_positive_area_gdf = gpd.GeoDataFrame(
        geometry=buffered_positive_area_gdf, crs=projected_crs
    )
    buffered_positive_area_gdf = buffered_positive_area_gdf[
        buffered_positive_area_gdf.geometry.notnull()
    ]
    buffered_positive_area_gdf = buffered_positive_area_gdf[
        buffered_positive_area_gdf.is_valid
    ]

    # Remove the overlapping areas (buffered) from the clipped shoreline
    clipped_shoreline_natura = gpd.overlay(
        clipped_shoreline, buffered_positive_area_gdf, how="difference"
    )

    # Ensure the final GeoDataFrames are in EPSG:4326 for plotting
    clipped_shoreline = clipped_shoreline.to_crs(epsg=4326)
    buffered_positive_area_gdf = buffered_positive_area_gdf.to_crs(epsg=4326)
    clipped_shoreline_natura = clipped_shoreline_natura.to_crs(epsg=4326)

    # Save to snakemake.output
    clipped_shoreline_natura.to_file(snakemake.output.shorelines_natura)
    buffered_positive_area_gdf.to_file(snakemake.output.buffered_natura)

    return clipped_shoreline_natura


def prepare_gebco(country_shapes, gebco_path):
    """
    Prepares GEBCO raster data by clipping it to country boundaries and reprojecting it.

    Parameters:
    - country_shapes (GeoDataFrame): GeoDataFrame containing country boundary geometries.
    - gebco_path (str): Path to the GEBCO raster file.

    Returns:
    - gebco_raster (rasterio.DatasetReader): Reprojected raster dataset clipped to the country boundaries.

    Functionality:
    - Clips the GEBCO raster using the unified geometry of the country shapes.
    - Saves the clipped raster to a file.
    - Reprojects the clipped raster to EPSG:3857 for distance calculations.
    - Saves the reprojected raster to a file.
    - Returns the reprojected raster for further processing.
    """

    # Clip the raster using the geometry of country_shapes
    with rasterio.open(gebco_path) as src:
        # Use the unified geometry of country_shapes
        country_geometry = country_shapes.geometry.unary_union

        # Clip the raster
        out_image, out_transform = mask(src, [country_geometry], crop=True)
        out_meta = src.meta.copy()
        out_meta.update(
            {
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform,
            }
        )

        # Save the clipped raster
        clipped_output_path = Path(
            os.path.join(BASE_DIR, "data/gebco/GEBCO_2024_clipped.tif")
        )
        with rasterio.open(clipped_output_path, "w", **out_meta) as dest:
            dest.write(out_image)

    # Reprojecct and save clipped raster
    reprojected_clipped_output_path = Path(
        os.path.join(BASE_DIR, "data/gebco/GEBCO_2024_reprojected_clipped.tif")
    )

    with rasterio.open(clipped_output_path) as src:
        # Target crs to be able to calculate distances
        target_crs = "EPSG:3857"

        # Calculate transformation and dimensions for the target CRS
        transform, width, height = calculate_default_transform(
            src.crs, target_crs, src.width, src.height, *src.bounds
        )

        # Update metadata
        kwargs = src.meta.copy()
        kwargs.update(
            {
                "crs": target_crs,
                "transform": transform,
                "width": width,
                "height": height,
            }
        )

        # Write the reprojected raster
        with rasterio.open(reprojected_clipped_output_path, "w", **kwargs) as dst:
            for i in range(1, src.count + 1):  # Iterate over raster bands
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=target_crs,
                    resampling=Resampling.nearest,
                )

    # Load the reprojected raster
    gebco_raster = rasterio.open(reprojected_clipped_output_path)

    return gebco_raster


def prepare_aqueduct(regions_onshore, aqueduct_file_path):
    """
    Prepares aqueduct data by filtering, intersecting, and aggregating water stress values.

    Parameters:
    - regions_onshore (GeoDataFrame): GeoDataFrame containing onshore regions with country information.
    - aqueduct_file_path (str): Path to the aqueduct GeoDatabase file.

    Returns:
    - regions_onshore_aqueduct (GeoDataFrame): GeoDataFrame with aggregated water stress values and additional columns for visualization.

    Functionality:
    - Filters aqueduct data to include only features matching the countries in `regions_onshore`.
    - Computes the intersection between aqueduct data and onshore regions.
    - Removes invalid water stress values and caps them to a maximum of 1.
    - Calculates the area of intersected geometries and aggregates water stress values using a weighted average.
    - Adds columns for scaled percentages, classification, and color mapping for visualization.
    - Saves the processed GeoDataFrame to a file specified by `snakemake.output.regions_onshore_aqueduct`.
    """
    # List all layers in the GeoDatabase
    # layers = gpd.io.file.fiona.listlayers(aqueduct_file_path)
    layer_name = "baseline_annual"  # TODO To put this in config file for user to choose if baseline or future

    # Create a column with iso3 country code
    cc = coco.CountryConverter()
    iso2 = pd.Series(regions_onshore.country)
    regions_onshore["iso3"] = cc.pandas_convert(series=iso2, to="ISO3")

    # Use a filter to read only features where gid_0 IN regions_onshore
    countries = regions_onshore.iso3.unique()
    country_list = "', '".join(countries)  # Join the country codes with commas
    filter_expression = f"gid_0 IN ('{country_list}')"

    # filter_expression"
    aqueduct = gpd.read_file(
        aqueduct_file_path, layer=layer_name, where=filter_expression
    )

    # Compute the intersection between aqueduct and regions_onshore
    intersected = gpd.overlay(regions_onshore, aqueduct, how="intersection")

    # Remove rows where `bws_raw` is outside [0, 1]
    intersected = intersected[(intersected["bws_raw"] >= 0)]
    # Cap values of `bws_raw` to the range [1]
    intersected["bws_raw"] = intersected["bws_raw"].clip(upper=1)

    # Calculate the area of each intersected geometry (handles MultiPolygons)
    intersected["area"] = intersected.to_crs(epsg="3857").geometry.area

    # Aggregate the `bws_raw` column using a weighted average based on area
    aggregated_values = (
        intersected.groupby("name")
        .apply(lambda x: (x["bws_raw"] * x["area"]).sum() / x["area"].sum())
        .reset_index(name="aggregated_bws_raw")
    )

    # Merge the aggregated values back into regions_onshore
    regions_onshore_aqueduct = regions_onshore.merge(
        aggregated_values, on="name", how="left"
    )

    # Scale values for plotting (convert to percentages)
    regions_onshore_aqueduct["aggregated_bws_raw_percent"] = (
        regions_onshore_aqueduct["aggregated_bws_raw"] * 100
    )

    # Create a new column to classify the data
    regions_onshore_aqueduct["color_category"] = regions_onshore_aqueduct[
        "aggregated_bws_raw_percent"
    ].apply(lambda x: "above_40" if x > 40 else "below_40")

    # Define a custom color map: Red for "above_40", Blue for "below_40"
    color_map = {"above_40": "red", "below_40": "blue"}
    regions_onshore_aqueduct["color"] = regions_onshore_aqueduct["color_category"].map(
        color_map
    )

    # Save to snakemake.output
    regions_onshore_aqueduct.to_file(snakemake.output.regions_onshore_aqueduct)

    return regions_onshore_aqueduct


def calc_distances_to_shore(
    regions_onshore_aqueduct, clipped_shoreline_natura, gebco_raster
):
    """
    Calculates distances and altitude changes from onshore aqueduct regions to the nearest shoreline.

    Parameters:
    - regions_onshore_aqueduct (GeoDataFrame): GeoDataFrame containing onshore aqueduct regions with water stress classification.
    - clipped_shoreline_natura (GeoDataFrame): GeoDataFrame containing shoreline geometries clipped by Natura 2000 areas.
    - gebco_raster (rasterio.DatasetReader): GEBCO raster dataset for altitude data.

    Returns:
    - shortest_distance_lines (GeoDataFrame): GeoDataFrame containing LineStrings connecting centroids to nearest shoreline points,
      along with calculated distances and altitude changes.

    Functionality:
    - Filters aqueduct regions with water stress classification 'above_40'.
    - Computes centroids for filtered regions and finds the nearest shoreline points.
    - Creates LineStrings connecting centroids to nearest shoreline points.
    - Calculates altitude differences along the LineStrings using GEBCO raster data.
    - Computes distances in kilometers and applies a routing factor for adjusted distances and altitude changes.
    - Reprojects geometries to EPSG:4326 for output.
    - Saves desalination regions to a file specified by `snakemake.output.regions_onshore_aqueduct_desalination`.
    """
    # Filter GeoDataFrame for regions with 'above_40'
    if snakemake.params.water_stress:
        gadm_desal = regions_onshore_aqueduct[
            regions_onshore_aqueduct.color_category == snakemake.params.water_stress
        ]

    else:
        gadm_desal = regions_onshore_aqueduct

    # Ensure both GeoDataFrames are in the same projected CRS
    if gadm_desal.crs.is_geographic:
        gadm_desal = gadm_desal.to_crs(epsg=3857)  # Example: Web Mercator (meters)
    if clipped_shoreline_natura.crs != gadm_desal.crs:
        clipped_shoreline_natura = clipped_shoreline_natura.to_crs(gadm_desal.crs)

    # Compute centroids and save them in a separate column
    gadm_desal["centroid"] = gadm_desal.geometry.centroid

    # Calculate the nearest point on the shoreline for each centroid
    def get_nearest_point(centroid, shoreline):
        # Find the closest geometry in the shoreline GeoDataFrame
        distances = shoreline.distance(centroid)
        nearest_geom = shoreline.iloc[distances.idxmin()].geometry

        # For MultiLineString or LineString, compute the true nearest point
        if nearest_geom.geom_type in ["MultiLineString", "LineString"]:
            # Use Shapely's `nearest_points` to find the exact nearest point
            nearest_point = nearest_points(centroid, nearest_geom)[1]
            return nearest_point
        else:
            raise ValueError(f"Unexpected geometry type: {nearest_geom.geom_type}")

    gadm_desal["nearest_point"] = gadm_desal["centroid"].apply(
        lambda centroid: get_nearest_point(centroid, clipped_shoreline_natura)
    )

    # Create LineStrings connecting centroids to their nearest points
    gadm_desal["shortest_distance_line"] = gadm_desal.apply(
        lambda row: LineString([row["nearest_point"], row["centroid"]]), axis=1
    )

    # # Create shortest_distance_lines GeoDataFrame with all columns from gadm_desal
    shortest_distance_lines = gadm_desal.drop(columns="geometry").copy()

    # Ensure the geometry of shortest_distance_lines is set to 'gadm_desal'
    shortest_distance_lines = shortest_distance_lines.set_geometry(
        "shortest_distance_line"
    )
    shortest_distance_lines = shortest_distance_lines.set_crs(gadm_desal.crs)

    # Extract bus names for centroids and nearest points
    def get_bus_name(point, regions):
        # Find the region where the point is within the geometry
        matching_region = regions[regions.geometry.intersects(point)]

        if not matching_region.empty:
            # Return the name of the first matching region
            return matching_region.iloc[0]["name"]
        else:
            raise ValueError("Point does not intersect with any region.")

    shortest_distance_lines["centroid_bus"] = (
        gadm_desal["centroid"]
        .to_crs(epsg=4326)
        .apply(lambda centroid: get_bus_name(centroid, regions_onshore_aqueduct))
    )

    # Create new GeoDataFrame using the original attributes, but set the new geometry column
    gdf_nearest = gadm_desal.copy()
    gdf_nearest.set_geometry("nearest_point", inplace=True)
    nearest = gpd.sjoin_nearest(
        gdf_nearest, gadm_desal, how="left", distance_col="dist"
    )

    shortest_distance_lines["nearest_point_bus"] = nearest["name_right"]

    shortest_distance_lines["line_name"] = (
        shortest_distance_lines["centroid_bus"]
        + " -> "
        + shortest_distance_lines["nearest_point_bus"]
    )

    # Apply the calculate_altitude_difference function to all shortest_distance_lines
    results = []
    dfs = []

    for line_name, line in zip(
        shortest_distance_lines.line_name, shortest_distance_lines.geometry
    ):
        df_i, total_altitude = calculate_altitude_difference(line, gebco_raster)
        df_i["line_name"] = line_name
        df_i = df_i.assign(**total_altitude)
        dfs.append(df_i)
        results.append(total_altitude["required_altitude"])

    shortest_distance_lines["total_positive_altitude_change"] = results
    all_profiles = pd.concat(dfs, ignore_index=True)

    # Compute distances and apply the routing factor
    shortest_distance_lines["distance_km"] = shortest_distance_lines.geometry.apply(
        calculate_distance_km
    )
    shortest_distance_lines["adjusted_distance_km"] = (
        shortest_distance_lines["distance_km"] * 1.5
    )
    shortest_distance_lines["adjusted_altitude_change"] = (
        shortest_distance_lines["total_positive_altitude_change"] / 1.5
    )

    # Ensure the columns are in the GeoDataFrame and reproject to EPSG:4326
    for col in ["centroid", "nearest_point"]:
        if (
            col in shortest_distance_lines.columns
            or shortest_distance_lines[col].dtype.name == "geometry"
        ):
            shortest_distance_lines[col] = shortest_distance_lines[col].to_crs(
                epsg=4326
            )
            shortest_distance_lines[col] = shortest_distance_lines[col].apply(
                lambda geom: geom.wkt if geom else None
            )

    # Ensure the columns are in the GeoDataFrame and reproject to EPSG:4326
    for col in ["centroid", "nearest_point", "shortest_distance_line"]:
        if col in gadm_desal.columns and gadm_desal[col].dtype.name == "geometry":
            gadm_desal[col] = GeoSeries(gadm_desal[col], crs=gadm_desal.crs).to_crs(
                epsg=4326
            )
            gadm_desal[col] = gadm_desal[col].apply(
                lambda geom: geom.wkt if geom else None
            )

    gadm_desal.to_file(snakemake.output.regions_onshore_aqueduct_desalination)

    return shortest_distance_lines, all_profiles


def calculate_altitude_difference(line, raster, resolution=1000, max_altitude=3000):
    """
    Sample elevation and compute hydrostatic running-head bookkeeping.

    Returns:
      df: DataFrame with columns:
          distance_m, elevation_m, diff_m (step), effective_increment_m,
          head_clamped_m (running), clipped_m (lost due to cap, per step)
      stats: dict with totals and start/end info
    """
    # --- sample distances (include endpoint) ---
    if line.length == 0:
        raise ValueError("Line has zero length.")
    d = np.arange(0, line.length + 1e-9, resolution)
    if d[-1] < line.length:
        d = np.append(d, line.length)

    pts = [line.interpolate(float(x)) for x in d]
    coords = [(p.x, p.y) for p in pts]

    # --- raster sampling (dataset method with fallback) ---
    try:
        it = raster.sample(coords)
    except AttributeError:
        it = sample.sample_gen(raster, coords)

    elev = np.array([float(v[0]) for v in it], dtype=float)

    # Enforce non-negative elevations (treat anything below sea level as 0 m)
    elev = np.clip(elev, 0.0, None)  # same as: elev[elev < 0] = 0.0

    # simple NaN handling (forward/back fill, then drop any remaining)
    if np.isnan(elev).any():
        for i in range(1, len(elev)):
            if np.isnan(elev[i]) and not np.isnan(elev[i - 1]):
                elev[i] = elev[i - 1]
        for i in range(len(elev) - 2, -1, -1):
            if np.isnan(elev[i]) and not np.isnan(elev[i + 1]):
                elev[i] = elev[i + 1]
        mask = ~np.isnan(elev)
        d, elev = d[mask], elev[mask]
        if len(elev) < 2:
            raise ValueError("Insufficient valid elevation samples after NaN handling.")

    # --- diffs and running-head bookkeeping ---
    diffs = np.diff(elev)
    head = 0.0
    head_series = [head]
    effective_inc = []
    clipped = []

    for diff in diffs:
        if diff > 0:
            room = max_altitude - head
            add = max(0.0, min(diff, room))  # what counts this step
            clip = max(0.0, diff - room)  # what is lost due to cap
            head = min(max_altitude, head + diff)
            effective_inc.append(add)
            clipped.append(clip)
        elif diff < 0:
            head = max(0.0, head + diff)
            effective_inc.append(0.0)
            clipped.append(0.0)
        else:
            effective_inc.append(0.0)
            clipped.append(0.0)
        head_series.append(head)

    head_series = np.array(head_series)
    effective_inc = np.array(effective_inc)
    clipped = np.array(clipped)

    elev_start, elev_end = float(elev[0]), float(elev[-1])
    min_required_head = max(elev_end - elev_start, 0.0)
    effective_head = float(effective_inc.sum())
    required_altitude = max(effective_head, min_required_head)

    df = pd.DataFrame(
        {"distance_m": d, "elevation_m": elev, "head_clamped_m": head_series}
    )
    # step-wise arrays align to the *end* of each segment
    step_df = pd.DataFrame(
        {
            "distance_m": d[1:],
            "diff_m": diffs,
            "effective_increment_m": effective_inc,
            "clipped_m": clipped,
        }
    )
    df = df.merge(step_df, on="distance_m", how="left")

    stats = {
        "effective_head": effective_head,
        "min_required_head": float(min_required_head),
        "required_altitude": float(required_altitude),
        "elev_start": elev_start,
        "elev_end": elev_end,
        "max_altitude_cap": float(max_altitude),
        "total_clipped": float(clipped.sum()),
        "resolution_m": float(resolution),
    }
    return df, stats


def calculate_distance_km(line):
    """
    Calculate the length of a LineString in kilometers.
    """
    return line.length / 1000  # Convert meters to kilometers


def delta_p_friction(lambda_f, L, d, rho, v):
    """
    Calculate pressure drop due to friction in a pipe.

    Parameters:
        lambda_f : float
            Friction factor [-]
        L : float
            Pipe length [m]
        d : float
            Pipe diameter [m]
        rho : float
            Fluid density [kg/m³]
        v : float
            Fluid velocity [m/s]

    Returns:
        delta_p_bar : float
            Pressure drop [bar]
    """
    delta_p_Pa = lambda_f * (L / d) * (rho * v**2 / 2)
    delta_p_bar = delta_p_Pa / 1e5
    return delta_p_bar


def delta_p_altitude(rho, h, g=9.81):
    """
    Calculate pressure drop due to altitude difference.

    Parameters:
        rho : float
            Fluid density [kg/m³]
        h : float
            Altitude difference [m]
        g : float, optional
            Gravitational acceleration [m/s²], default is 9.81.

    Returns:
        delta_p_Pa : float
            Pressure difference [Pa]
        delta_p_bar : float
            Pressure difference [bar]
    """
    delta_p_Pa = rho * g * h
    delta_p_bar = delta_p_Pa / 1e5
    return delta_p_bar


def number_of_pumping_stations(delta_p_total, p_max_pipeline):
    """
    Calculate the required number of pumping stations.

    Parameters:
        delta_p_total : float
            Total pressure drop [bar]
        p_max_pipeline : float
            Nominal pressure of the pipeline [bar]

    Returns:
        n : float
            Number of pumping stations needed
    """
    n = delta_p_total / p_max_pipeline
    # return round(n,0)
    return n


def pump_electrical_power_kw(
    total_dp, Nbr_of_pumping_stations, mass_flow_rate, rho, eta_pump, eta_motor
):
    """
    Calculate the electrical power consumption of a pump in kW.

    Parameters:
        total_dp : float
            Total pressure drop [bar]
        Nbr_of_pumping_stations: int
            Number of pumping stations [-]
        mass_flow_rate : float
            Mass flow rate [kg/s]
        rho : float
            Fluid density [kg/m³]
        eta_pump : float
            Pump efficiency [-], e.g., 0.85
        eta_motor : float
            Motor efficiency [-], e.g., 0.92

    Returns:
        P_el_pump_W : float
            Electrical power [W]
        P_el_pump_kW : float
            Electrical power [kW]
    """
    delta_p_pump = (total_dp / Nbr_of_pumping_stations) * 1e5  # Convert bar to Pa
    P_el_pump_W = delta_p_pump * (mass_flow_rate / rho) * (1 / (eta_pump * eta_motor))
    P_el_pump_kW = P_el_pump_W / 1000  # Convert W → kW
    # return round(P_el_pump_kW,1)
    return P_el_pump_kW


def pump_station_cost_kw(power_kw, usd_to_eur, a=35768, b=0.558):
    """
    Calculate pump station investment cost from installed power in kW.

    Parameters:
    -----------
    power_kw : float or array
        Installed pump power in kW.
    a : float
        Cost factor in Eur/HP (default 35768). Source: Inflation a2023=15570(1.02)**42≈35,768 USD/HP in https://hypat.de/hypat-wAssets/docs/new/publikationen/HYPAT_WP_Water-Supply-for-Electrolysis-Plants.pdf).
    b : float
        Scaling exponent (default 0.558) Source: 1+b = 1-0.442 in https://hypat.de/hypat-wAssets/docs/new/publikationen/HYPAT_WP_Water-Supply-for-Electrolysis-Plants.pdf).
    usd_to_eur : float
        Conversion rate from USD to EUR .

    Returns:
    --------
    cost_eur : float or array
        Pump station cost in million EUR.
    """
    # Convert kW to HP
    power_hp = power_kw / 0.7457

    # Cost formula in USD
    cost_usd = a * (power_hp**b)

    # Convert to EUR
    cost_eur = cost_usd * usd_to_eur

    # Return cost in million EUR
    return cost_eur / 1e6


def add_pipeline_hydraulics(row):
    """
    Calculate pipeline hydraulics, pumping station requirements, and associated costs.

    This function computes various hydraulic parameters, including pressure drops due to friction
    and altitude, the number of pumping stations required, pump electrical power, and investment
    costs for both pumps and pipelines. It also calculates electricity demand and efficiency metrics.

    Parameters:
    -----------
    row : pd.Series
        A pandas Series containing the following keys:
        - 'adjusted_distance_km': float, the adjusted pipeline distance in kilometers.
        - 'adjusted_altitude_change': float, the altitude difference in meters.

    Returns:
    --------
    pd.Series
        A pandas Series containing the following calculated values:
        - 'dp_friction_bar': float, pressure drop due to friction in bar.
        - 'dp_altitude_bar': float, pressure drop due to altitude in bar.
        - 'total_dp_bar': float, total pressure drop in bar.
        - 'n_pumping_stations': int, number of pumping stations required.
        - 'power_kW': float, pump electrical power in kW.
        - 'mass_flow_rate_m3h': float, mass flow rate in m³/h.
        - 'invest_pumping_station': float, investment cost per pumping station in million EUR.
        - 'total_pump_invest_mil_eur': float, total pump investment cost in million EUR.
        - 'total_pump_invest_eur': float, total pump investment cost in EUR.
        - 'total_pump_invest_eur_per_kW': float, pump investment cost per kW in EUR.
        - 'pump_elec_demand': float, pump electricity demand in MWh per m³.
        - 'pump_capex_per_pnom': float, pump CAPEX per unit flow (€/m³/h).
        - 'pipeline_capex_per_pnom': float, pipeline CAPEX per unit flow (€/m³/h).
        - 'efficiency2': float, efficiency metric (-kWh_el per m³).

    Notes:
    ------
    - The function assumes a constant mass flow rate of 34.71 kg/s, corresponding to a 1000 MW electrolyzer.
    - The pipeline investment formulas are based on HYPAT formulas.
    - The pump station cost is adjusted for inflation from 1981 to 2023.

    References:
    -----------
    - HYPAT formulas: https://hypat.de/hypat-wAssets/docs/new/publikationen/HYPAT_WP_Water-Supply-for-Electrolysis-Plants.pdf
    """

    # Constansts
    d_mm = 171.65
    d_m = d_mm / 1000
    v = 1.5
    rho = 1000
    eta_pump, eta_motor = 0.85, 0.90
    eta_sys = eta_pump * eta_motor
    lambda_f = 0.01489
    mass_flow_rate = 34.71  # kg/s, Considered mass flow rate for 1000 MW Electrolyzer
    mass_flow_rate_m3s = mass_flow_rate / rho
    mass_flow_rate_m3h = mass_flow_rate_m3s * 3600

    # Calculate pressure drop due to friction
    dp_friction = delta_p_friction(
        lambda_f=lambda_f,  # example friction factor
        L=row["adjusted_distance_km"] * 1000,  # length in meters
        d=d_m,  # diameter in meters
        rho=rho,  # density in kg/m³
        v=v,  # velocity in m/s
    )

    # Calculate pressure drop due to altitude
    dp_altitude = delta_p_altitude(
        rho=rho,  # water density in kg/m³
        h=row["adjusted_altitude_change"],  # altitude difference in meters
    )

    # Total pressure drop
    total_dp = dp_friction + dp_altitude

    # Number of pumping stations
    n_pumping_stations = number_of_pumping_stations(
        delta_p_total=total_dp,
        p_max_pipeline=16,  # nominal pressure of the pipeline in bar
    )

    # Pump electrical power
    power_kW = pump_electrical_power_kw(
        total_dp=total_dp,
        Nbr_of_pumping_stations=n_pumping_stations,
        mass_flow_rate=mass_flow_rate,  # kg/s
        rho=rho,  # kg/m³
        eta_pump=eta_pump,  # pump efficiency
        eta_motor=eta_motor,  # motor efficiency
    )

    # Pump investment in millions Eur
    invest_pumping_station = pump_station_cost_kw(
        power_kw=power_kW,
        usd_to_eur=snakemake.params.costs["default_exchange_rate"],
        a=35768,  # 15570 # 35768 with inflation 2% 1981 till 2023
        b=0.558,
    )

    # Pipeline investment (HYPAT formulas: https://hypat.de/hypat-wAssets/docs/new/publikationen/HYPAT_WP_Water-Supply-for-Electrolysis-Plants.pdf)
    C_pipes = 1.1852 * (d_mm**1.9557)  # €/km material
    f_inst = 176.97 * (d_mm**-0.624)  # installation factor
    C_inst_km = C_pipes * f_inst  # €/km installed
    # Pipeline investment per km per m3/h
    capex_per_pnom_per_km = C_inst_km / mass_flow_rate_m3h  # €/ (m3/h) / km
    L_km = row["adjusted_distance_km"]  # length in km

    # pipeline CAPEX per unit flow (€/ (m3/h))
    pipeline_capex_per_pnom = capex_per_pnom_per_km * L_km

    # pump CAPEX per unit flow (€/ (m3/h))
    total_pump_invest_eur = invest_pumping_station * n_pumping_stations * 1e6
    pump_capex_per_pnom = total_pump_invest_eur / mass_flow_rate_m3h

    # --- electricity per m3 ---
    dp_fric_Pa = dp_friction * 1e5
    e_fric = dp_fric_Pa / (eta_sys * 3.6e6)  # kWh/m3
    dp_alt_Pa = dp_altitude * 1e5
    e_alt = max(0.0, dp_alt_Pa / (eta_sys * 3.6e6))

    efficiency2 = -(e_fric + e_alt) / 1000  # MW per (m³/h)

    pump_elec_demand = (e_fric + e_alt) / 1000  # MWh per m³

    return pd.Series(
        {
            "dp_friction_bar": dp_friction,
            "dp_altitude_bar": dp_altitude,
            "total_dp_bar": total_dp,
            "n_pumping_stations": n_pumping_stations,
            "power_kW": power_kW,
            "mass_flow_rate_m3h": mass_flow_rate_m3h,
            "invest_pumping_station": invest_pumping_station,
            "total_pump_invest_mil_eur": invest_pumping_station * n_pumping_stations,
            "total_pump_invest_eur": invest_pumping_station * n_pumping_stations * 1e6,
            "total_pump_invest_eur_per_kW": invest_pumping_station
            * 1e6
            / power_kW,  # Convert to Eur/kW
            "pump_elec_demand": pump_elec_demand,  # MWh_el per m3
            "pump_capex_per_pnom": pump_capex_per_pnom,
            "pipeline_capex_per_pnom": pipeline_capex_per_pnom,
            "efficiency2": efficiency2,  # -kWh_el per m3
        }
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_water_network",
            simpl="",
            clusters="30",
        )

    shoreline_gdf = download_shorelines().copy()
    gebco_path = download_gebco()
    aqueduct_file_path = download_aqueduct()

    # Load input data
    country_shapes = gpd.read_file(snakemake.input.country_shapes)
    regions_onshore = gpd.read_file(snakemake.input.regions_onshore)
    natura_tiff_path = snakemake.input.natura

    # Clip shoreline to country borders
    clipped_shoreline = clip_shorelines_country(shoreline_gdf, country_shapes)

    # Exclude Natura
    clipped_shoreline_natura = clip_shorelines_natura(
        natura_tiff_path, country_shapes, clipped_shoreline
    )

    # Prepare Digital Elevation Raster
    gebco_raster = prepare_gebco(country_shapes, gebco_path)

    # Prepare Water Stress regions (Aqueduct)
    regions_onshore_aqueduct = prepare_aqueduct(regions_onshore, aqueduct_file_path)

    # Calculate Distances and Elevations to shore
    shortest_distance_lines, water_pipes_profiles = calc_distances_to_shore(
        regions_onshore_aqueduct, clipped_shoreline_natura, gebco_raster
    )

    # Add pipeline hydraulics calculations
    hydraulics = shortest_distance_lines.apply(add_pipeline_hydraulics, axis=1)

    # Combine the results with the original GeoDataFrame
    shortest_distance_lines = shortest_distance_lines.join(hydraulics)

    water_pipes_profiles.to_csv(snakemake.output.water_pipes_profiles)

    # Save to snakemake.output
    shortest_distance_lines.to_file(snakemake.output.clustered_water_network)
