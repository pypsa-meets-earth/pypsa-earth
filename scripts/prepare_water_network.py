# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Prepare water network for Desalination and input to Electrolysis in prepare_sector_network.
"""

import logging

logger = logging.getLogger(__name__)

import os
import geopandas as gpd
import pandas as pd
import requests
import zipfile
import io
import matplotlib.pyplot as plt
import rasterio
from rasterio.plot import show
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from matplotlib.colors import LinearSegmentedColormap, to_rgba
import matplotlib.patches as mpatches
from shapely.geometry import LineString
from rasterio import sample
import country_converter as coco
import numpy as np
from rasterio.features import shapes
from shapely.geometry import shape
from shapely.ops import nearest_points
from geopandas import GeoSeries
from pathlib import Path
import xarray as xr
from _helpers import (
    BASE_DIR,
    content_retrieve,
    progress_retrieve,
    two_2_three_digits_country,
)


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
        logger.info(f"Gebco data already exists at '{gebco_file_path}'. Skipping download.")
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
        logger.info(f"Aqueduct data already exists at '{aqueduct_file_path}'. Skipping download.")
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
    url = "https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-2.3.7.zip"
    shoreline_file = "GSHHS_shp/i/GSHHS_i_L1.shp"  # Relative path inside the zip
    output_file = snakemake.output.shorelines

    # Save locations
    zip_fn = Path(os.path.join(BASE_DIR, "shorelines.zip"))
    to_fn = Path(os.path.join(BASE_DIR, "data/shorelines"))

    # Path to the extracted shapefile
    extracted_shoreline_path = to_fn / shoreline_file

    # Check if the target file already exists
    if extracted_shoreline_path.exists():
        logger.info(f"Shoreline data already exists at '{extracted_shoreline_path}'. Skipping download.")

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
    country_gdf['geometry'] = country_gdf.buffer(3000)

    # Convert back to the original CRS if needed
    country_gdf = country_gdf.to_crs(epsg=4326)

    # Ensure CRS consistency
    if shoreline_gdf.crs != country_gdf.crs:
        shoreline_gdf = shoreline_gdf.to_crs(country_gdf.crs)

    # Extract the boundary (exterior) of the shoreline polygons
    shoreline_boundary = shoreline_gdf.boundary

    # Convert the boundary to a GeoDataFrame
    shoreline_boundary_gdf = gpd.GeoDataFrame(geometry=shoreline_boundary, crs=shoreline_gdf.crs)

    # Clip the shoreline boundary to the country boundary
    clipped_shoreline = gpd.overlay(shoreline_boundary_gdf, country_gdf, how="intersection")

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
    positive_area_gdf = gpd.GeoDataFrame({'geometry': positive_polygons}, crs=src.crs)

    # Switch to a projected CRS for buffering (e.g., EPSG:3857)
    projected_crs = "EPSG:3857"
    positive_area_gdf = positive_area_gdf.to_crs(projected_crs)
    clipped_shoreline = clipped_shoreline.to_crs(projected_crs)

    # Buffer the positive areas by 4 km (4000 meters)
    buffered_positive_area_gdf = positive_area_gdf.buffer(4000)

    # Drop empty or invalid geometries after buffering
    buffered_positive_area_gdf = gpd.GeoDataFrame(geometry=buffered_positive_area_gdf, crs=projected_crs)
    buffered_positive_area_gdf = buffered_positive_area_gdf[buffered_positive_area_gdf.geometry.notnull()]
    buffered_positive_area_gdf = buffered_positive_area_gdf[buffered_positive_area_gdf.is_valid]

    # Remove the overlapping areas (buffered) from the clipped shoreline
    clipped_shoreline_natura = gpd.overlay(clipped_shoreline, buffered_positive_area_gdf, how="difference")

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
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform
        })

        # Save the clipped raster
        clipped_output_path = Path(os.path.join(BASE_DIR, "data/gebco/GEBCO_2024_clipped.tif"))
        with rasterio.open(clipped_output_path, "w", **out_meta) as dest:
            dest.write(out_image)
        
    # Reprojecct and save clipped raster
    reprojected_clipped_output_path = Path(os.path.join(BASE_DIR, "data/gebco/GEBCO_2024_reprojected_clipped.tif"))

    with rasterio.open(clipped_output_path) as src:
        # Target crs to be able to calculate distances
        target_crs = "EPSG:3857" 

        # Calculate transformation and dimensions for the target CRS
        transform, width, height = calculate_default_transform(
            src.crs, target_crs, src.width, src.height, *src.bounds
        )

        # Update metadata
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': target_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        # Write the reprojected raster
        with rasterio.open(reprojected_clipped_output_path, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):  # Iterate over raster bands
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=target_crs,
                    resampling=Resampling.nearest
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
    layers = gpd.io.file.fiona.listlayers(aqueduct_file_path)
    layer_name = 'baseline_annual'  # TODO To put this in config file for user to choose if baseline or future

    # Create a column with iso3 country code
    cc = coco.CountryConverter()
    iso2 = pd.Series(regions_onshore.country)
    regions_onshore["iso3"] = cc.pandas_convert(series=iso2, to='ISO3') 

    # Use a filter to read only features where gid_0 IN regions_onshore
    countries = regions_onshore.iso3.unique()
    country_list = "', '".join(countries)  # Join the country codes with commas
    filter_expression = f"gid_0 IN ('{country_list}')"

    # filter_expression"
    aqueduct = gpd.read_file(aqueduct_file_path, layer=layer_name, where=filter_expression)


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
    regions_onshore_aqueduct = regions_onshore.merge(aggregated_values, on="name", how="left")

    # Scale values for plotting (convert to percentages)
    regions_onshore_aqueduct['aggregated_bws_raw_percent'] = regions_onshore_aqueduct['aggregated_bws_raw'] * 100

    # Create a new column to classify the data
    regions_onshore_aqueduct['color_category'] = regions_onshore_aqueduct['aggregated_bws_raw_percent'].apply(
        lambda x: 'above_40' if x > 40 else 'below_40'
    )

    # Define a custom color map: Red for "above_40", Blue for "below_40"
    color_map = {'above_40': 'red', 'below_40': 'blue'}
    regions_onshore_aqueduct['color'] = regions_onshore_aqueduct['color_category'].map(color_map)

    # Save to snakemake.output
    regions_onshore_aqueduct.to_file(snakemake.output.regions_onshore_aqueduct)

    return regions_onshore_aqueduct


def calc_distances_to_shore (regions_onshore_aqueduct, clipped_shoreline_natura, gebco_raster):
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
    gadm_desal = regions_onshore_aqueduct[regions_onshore_aqueduct.color_category == 'above_40']

    # Ensure both GeoDataFrames are in the same projected CRS
    if gadm_desal.crs.is_geographic:
        gadm_desal = gadm_desal.to_crs(epsg=3857)  # Example: Web Mercator (meters)
    if clipped_shoreline_natura.crs != gadm_desal.crs:
        clipped_shoreline_natura = clipped_shoreline_natura.to_crs(gadm_desal.crs)

    # Compute centroids and save them in a separate column
    gadm_desal['centroid'] = gadm_desal.geometry.centroid

    # Calculate the nearest point on the shoreline for each centroid
    def get_nearest_point(centroid, shoreline):
        # Find the closest geometry in the shoreline GeoDataFrame
        distances = shoreline.distance(centroid)
        nearest_geom = shoreline.iloc[distances.idxmin()].geometry

        # For MultiLineString or LineString, compute the true nearest point
        if nearest_geom.geom_type in ['MultiLineString', 'LineString']:
            # Use Shapely's `nearest_points` to find the exact nearest point
            nearest_point = nearest_points(centroid, nearest_geom)[1]
            return nearest_point
        else:
            raise ValueError(f"Unexpected geometry type: {nearest_geom.geom_type}")

    gadm_desal['nearest_point'] = gadm_desal['centroid'].apply(
        lambda centroid: get_nearest_point(centroid, clipped_shoreline_natura)
    )

    # Create LineStrings connecting centroids to their nearest points
    gadm_desal['shortest_distance_line'] = gadm_desal.apply(
        lambda row: LineString([row['centroid'], row['nearest_point']]), axis=1
    )


    
    # # Create shortest_distance_lines GeoDataFrame with all columns from gadm_desal
    shortest_distance_lines = gadm_desal.drop(columns='geometry').copy()

    # Ensure the geometry of shortest_distance_lines is set to 'gadm_desal'
    shortest_distance_lines = shortest_distance_lines.set_geometry('shortest_distance_line')
    shortest_distance_lines = shortest_distance_lines.set_crs(gadm_desal.crs)

    # # Calculate the shortest distance from centroids to the clipped shoreline
    # gadm_desal['shortest_distance_m'] = gadm_desal['centroid'].apply(
    #     lambda point: clipped_shoreline_natura.distance(point).min()
    # )

    # # Convert the distance from meters to kilometers
    # gadm_desal['shortest_distance_km'] = gadm_desal['shortest_distance_m'] / 1000


    # Apply the calculate_altitude_difference function to all shortest_distance_lines
    results = []

    for line in shortest_distance_lines.geometry:
        total_altitude = calculate_altitude_difference(line, gebco_raster)
        results.append(total_altitude)

    # Add the results to the GeoDataFrame
    shortest_distance_lines['total_positive_altitude_change'] = results

    # Compute distances and apply the routing factor
    shortest_distance_lines['distance_km'] = shortest_distance_lines.geometry.apply(calculate_distance_km)
    shortest_distance_lines['adjusted_distance_km'] = shortest_distance_lines['distance_km'] * 1.5
    shortest_distance_lines['adjusted_altitude_change'] = shortest_distance_lines['total_positive_altitude_change'] / 1.5

    # Ensure the columns are in the GeoDataFrame and reproject to EPSG:4326
    for col in ['centroid', 'nearest_point']:
        if col in shortest_distance_lines.columns or shortest_distance_lines[col].dtype.name == 'geometry':
            shortest_distance_lines[col] = shortest_distance_lines[col].to_crs(epsg=4326)
            shortest_distance_lines[col] = shortest_distance_lines[col].apply(lambda geom: geom.wkt if geom else None)

    # Ensure the columns are in the GeoDataFrame and reproject to EPSG:4326
    for col in ['centroid', 'nearest_point', 'shortest_distance_line']:
        if col in gadm_desal.columns and gadm_desal[col].dtype.name == 'geometry':
            gadm_desal[col] = GeoSeries(gadm_desal[col], crs=gadm_desal.crs).to_crs(epsg=4326)
            gadm_desal[col] = gadm_desal[col].apply(lambda geom: geom.wkt if geom else None)

    gadm_desal.to_file(snakemake.output.regions_onshore_aqueduct_desalination)


    return shortest_distance_lines






def calculate_altitude_difference(line, raster, resolution=1000, max_altitude=160):
    """
    Calculate the total positive change in altitude along a line, considering
    hydrostatic pressure constraints.

    Parameters:
    - line: Shapely LineString object.
    - raster: Rasterio dataset object.
    - resolution: Sampling resolution in meters (default: 1 km).
    - max_altitude: Maximum altitude difference limit (default: 160 m).

    Returns:
    - total_positive_altitude_change: Total positive altitude difference along the line.
    """
    # Sample points along the line every `resolution` meters
    distances = np.arange(0, line.length, resolution)
    sampled_points = [line.interpolate(d) for d in distances]

    # Extract raster values at sampled points
    coords = [(point.x, point.y) for point in sampled_points]
    raster_values = [
        val[0] for val in sample.sample_gen(raster, coords) if val[0] is not None
    ]

    # Calculate altitude differences between consecutive points
    altitude_differences = np.diff(raster_values)

    # Process altitude differences considering constraints
    total_positive_altitude_change = 0
    buffer = 0  # Tracks the impact of negative differences on the next segment

    for diff in altitude_differences:
        if diff > 0:  # Positive altitude difference
            # Apply buffer from previous negative differences
            diff = max(0, diff + buffer)
            # Apply the maximum altitude constraint
            diff = min(diff, max_altitude)
            total_positive_altitude_change += diff
            buffer = 0  # Reset buffer
        elif diff < 0:  # Negative altitude difference
            buffer += diff  # Accumulate buffer (negative values reduce positive altitude)

    return total_positive_altitude_change



def calculate_distance_km(line):
    """
    Calculate the length of a LineString in kilometers.
    """
    return line.length / 1000  # Convert meters to kilometers





if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_water_network",
            simpl="",
            clusters="10",
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
    clipped_shoreline_natura = clip_shorelines_natura(natura_tiff_path, country_shapes, clipped_shoreline)

    # Prepare Digital Elevation Raster
    gebco_raster = prepare_gebco(country_shapes, gebco_path)

    # Prepare Water Stress regions (Aqueduct)
    regions_onshore_aqueduct = prepare_aqueduct(regions_onshore, aqueduct_file_path)

    # Calculate Distances and Elevations to shore
    shortest_distance_lines = calc_distances_to_shore (regions_onshore_aqueduct, clipped_shoreline_natura, gebco_raster)
    
    # Save to snakemake.output
    shortest_distance_lines.to_file(snakemake.output.clustered_water_network)

