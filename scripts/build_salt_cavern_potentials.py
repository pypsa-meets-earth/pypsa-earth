# -*- coding: utf-8 -*-
"""
Build salt cavern potentials for hydrogen storage.

https://doi.org/10.3133/sir20105090S
"""

import functools
import math
import os
import shutil
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio.features
import requests
import rioxarray
import shapely.geometry
from _helpers import (
    mock_snakemake,
    to_csv_nafix,
    COPERNICUS_CRS
)

# TODO: consider externalizing this download in a retrieve_* rule or helper function
def download_potash_data():
    # URL of the Potash GIS data
    url = "https://pubs.usgs.gov/sir/2010/5090/s/PotashGIS.zip"

    # Download directory
    download_dir = "data/potash_gis"
    zip_path = os.path.join(download_dir, "PotashGIS.zip")

    # Ensure directory exists
    os.makedirs(download_dir, exist_ok=True)

    # Download the ZIP
    response = requests.get(url)
    if response.status_code == 200:
        with open(zip_path, "wb") as f:
            f.write(response.content)
    else:
        raise Exception(f"Download failed with status code {response.status_code}")

    # Extract full archive to preserve directory structure
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(path=download_dir)

    # Remove ZIP after extraction
    os.remove(zip_path)

    # Search for the Shapefile (.shp) within the extracted folder
    shp_path = None
    for root, dirs, files in os.walk(download_dir):
        for file in files:
            if file == "PotashTracts.shp":
                shp_path = os.path.join(root, file)
                break
        if shp_path:
            break

    if not shp_path:
        raise FileNotFoundError("PotashTracts.shp not found after extraction.")

    return gpd.read_file(shp_path)


def capsule_volume(diameter_m, height_m):
    """
    Calculate the volume of a cylindrical capsule-shaped cavern.

    A typical salt cavern for hydrogen storage is approximated as a capsule:
    a cylinder with hemispherical ends. This function computes the total
    volume based on the given diameter and total height (including the domed ends).

    Parameters
    ----------
    diameter_m : float
        Diameter of the cavern in meters.
    height_m : float
        Total height of the cavern in meters, including the hemispherical ends.

    Returns
    -------
    float
        Total cavern volume in cubic meters (m³).

    Notes
    -----
    The formula used is:
        V = π * (d / 2)² * (h - d) + (4/3) * π * (d / 2)³

    where:
        d = diameter
        h = total height
        (h - d) = cylindrical height
        two hemispheres together form a full sphere of diameter d
    """
    radius = diameter_m / 2
    h_cyl = height_m - diameter_m
    volume_cylinder = math.pi * radius**2 * h_cyl
    volume_spheres = (4 / 3) * math.pi * radius**3
    return volume_cylinder + volume_spheres


def compute_physical_capacity(
    diameter_m, height_m, rho_min, rho_max, theta_safety, lhv_h2
):
    """
    Computes the physical hydrogen storage energy capacity of an underground cavern in GWh.

    Parameters
    ----------
    diameter_m : float
        Diameter of the cavern in meters.
    height_m : float
        Height of the cavern in meters.
    rho_min : float
        Minimum hydrogen density in the cavern (e.g., at minimum pressure) in kg/m³.
    rho_max : float
        Maximum hydrogen density in the cavern (e.g., at maximum pressure) in kg/m³.
    theta_safety : float
        Usable fraction of the working gas capacity (typically < 1 to account for safety margins).
    lhv_h2 : float
        Lower heating value (LHV) of hydrogen in MJ/kg.

    Returns
    -------
    float
        Physical energy capacity of the cavern in GWh.

    Notes
    -----
    This function assumes a capsule-shaped cavern geometry and uses the LHV of hydrogen to
    compute the usable energy content of the working gas volume.
    """
    v_cavern = capsule_volume(diameter_m, height_m)
    m_working = (rho_max - rho_min) * v_cavern * theta_safety
    return m_working * lhv_h2 / 1e6  # GWh


def compute_gwh_per_km2(
    diameter_m,
    height_m,
    energy_density_MWh_per_m3=0.045,
    surface_area_per_cavern_km2=13.0,
):
    """
    Compute the hydrogen storage potential in GWh per km² for a capsule-shaped salt cavern.

    Parameters
    ----------
    diameter_m
        Diameter of the cavern in meters.
    height_m
        Total height of the cavern in meters.
    energy_density_MWh_per_m3, optional
        Usable energy content per cubic meter of hydrogen (default is 0.045 MWh/m³).
    surface_area_per_cavern_km2, optional
        Surface area required per cavern in km² (default is 13 km²).

    Returns
    -------
    Hydrogen storage potential in GWh per km² of surface area.
    """
    volume_m3 = capsule_volume(diameter_m, height_m)
    energy_per_cavern_MWh = volume_m3 * energy_density_MWh_per_m3
    caverns_per_km2 = 1 / surface_area_per_cavern_km2
    return (energy_per_cavern_MWh * caverns_per_km2) / 1000  # Convert to GWh


def classify_salt_type(gdf):
    """
    Classify the type of salt deposit based on geological and deposit type information.

    This function adds a new column 'salt_type' to the input GeoDataFrame, classifying
    each entry into one of the following categories:

    - "excluded_brine": if the 'Dep_type' column contains the word "brine"
    - "dome": if the 'Dep_type' or 'Geology' column contains the word "halokinetic"
    - "bedded": if the 'Dep_type' or 'Geology' column contains the word "stratabound" or "evaporite"
    - "unknown": if none of the above conditions are met

    Parameters:
        gdf (GeoDataFrame): A GeoPandas GeoDataFrame containing at least the columns
            'Dep_type' and 'Geology'.

    Returns:
        GeoDataFrame: The input GeoDataFrame with an additional column 'salt_type'
        indicating the classified salt deposit type.
    """

    def classify(row):
        dep = str(row["Dep_type"]).lower()
        geo = str(row["Geology"]).lower() if pd.notnull(row["Geology"]) else ""

        if "brine" in dep:
            return "excluded_brine"
        if "halokinetic" in dep or "halokinetic" in geo:
            return "dome"
        if (
            "stratabound" in dep
            or "evaporite" in dep
            or "stratabound" in geo
            or "evaporite" in geo
        ):
            return "bedded"
        return "unknown"

    gdf["salt_type"] = gdf.apply(classify, axis=1)
    return gdf["salt_type"]


def apply_landuse_exclusions(gdf, regions):
    """
    Applies land use exclusion zones to a GeoDataFrame of underground salt storage areas
    by removing regions intersecting with specified Copernicus land use types.

    The function performs the following steps:
    1. Applies a negative buffer to each salt storage geometry, with buffer size depending
       on whether the salt type is 'bedded' or 'dome'.
    2. Loads Copernicus land cover data and extracts raster cells matching specified grid codes.
    3. Converts matching raster regions into vector geometries and applies optional buffers.
    4. Merges all exclusion geometries and removes intersecting salt storage areas.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing underground salt storage areas with a column 'salt_type'
        and geometries in any CRS.
    area_crs : str or pyproj.CRS
        Target coordinate reference system (CRS) in which geometries will be processed and returned.

    Returns
    -------
    geopandas.GeoDataFrame
        A filtered GeoDataFrame in the specified CRS, excluding areas that intersect with
        Copernicus land use exclusions.

    Notes
    -----
    - Buffers for 'bedded' and 'dome' salt types, as well as Copernicus grid codes and
      optional per-code buffer distances, are loaded from `snakemake.params.underground_storage`.
    - This function assumes that `paths.copernicus` and `distance_crs` are defined elsewhere
      in the Snakemake workflow context.
    """
    cop = snakemake.params.underground_storage["copernicus"]
    salt_buffer_m_bedded = snakemake.params.underground_storage["salt_buffer_m_bedded"]
    salt_buffer_m_dome = snakemake.params.underground_storage["salt_buffer_m_dome"]

    gdf = gdf.to_crs(distance_crs).copy()

    # Apply buffer by salt type before masking
    gdf["geometry"] = gdf.apply(
        lambda row: row.geometry.buffer(
            -salt_buffer_m_bedded if row.salt_type == "bedded" else -salt_buffer_m_dome
        ),
        axis=1,
    )
    gdf = gdf[gdf.is_valid & ~gdf.is_empty]

    exclusion_shapes = []

    # Loading Copernicus raster
    with rioxarray.open_rasterio(paths.copernicus, masked=True) as rds:
        rds_clip = rds.squeeze().rio.clip_box(*regions.total_bounds)
        rds_reproj = rds_clip.rio.reproject(distance_crs)
        clc = rds_reproj.load()
        transform = clc.rio.transform()

        for code in cop["grid_codes"]:
            # Build mask
            mask = clc.data == code
            if not np.any(mask):
                continue

            shapes_gen = rasterio.features.shapes(
                mask.astype(np.uint8), mask=mask, transform=transform
            )

            # Extract geometries where value == 1
            code_shapes = [
                shapely.geometry.shape(s) for s, val in shapes_gen if val == 1
            ]

            if not code_shapes:
                continue

            gdf_mask = gpd.GeoDataFrame(geometry=code_shapes, crs=distance_crs)

            # Apply buffer if specified
            buffer_dist = (cop.get("distance_buffer_by_code") or {}).get(str(code), 0)
            if buffer_dist > 0:
                gdf_mask["geometry"] = gdf_mask.buffer(buffer_dist)
                gdf_mask = gdf_mask[gdf_mask.is_valid & ~gdf_mask.is_empty]

            exclusion_shapes.append(gdf_mask)

    if not exclusion_shapes:
        # No Copernicus exclusion zones found. Returning original GeoDataFrame.
        return gdf.to_crs(geo_crs)

    # Merging all exclusion zones
    all_exclusions = pd.concat(exclusion_shapes, ignore_index=True)
    union_geom = all_exclusions.geometry.unary_union
    exclusion_union = gpd.GeoDataFrame(geometry=[union_geom], crs=distance_crs)

    # Removing excluded regions
    gdf["geometry"] = gdf.geometry.difference(exclusion_union.geometry.iloc[0])
    gdf = gdf[gdf.is_valid & ~gdf.is_empty]

    return gdf.to_crs(geo_crs)


def estimate_h2_potential_from_potash(potash_gdf, regions, min_area_km2=13.0):
    """
    Estimate technical underground hydrogen storage potential from potash tracts,
    based on geological characteristics and cavern design assumptions.

    Parameters
    ----------
    potash_gdf : GeoDataFrame
        Raw potash GIS data.
    min_area_km2 : float
        Minimum area for a salt structure to be considered.

    Returns
    -------
    GeoDataFrame
        Filtered and annotated GeoDataFrame with salt cavern potential (GWh).
    """
    gdf = potash_gdf.copy()
    gdf["salt_type"] = classify_salt_type(gdf)

    # filter for valid salt types and minimum area
    valid_types = ["bedded", "dome"]
    filtered = gdf[gdf["salt_type"].isin(valid_types)].copy()

    filtered = apply_landuse_exclusions(filtered, regions)

    # Filter after exclusions to remove tiny remnants
    filtered = filtered[filtered["Area_km2"] >= min_area_km2]

    method = snakemake.params.underground_storage["method"]

    energy_density_MWh_per_m3 = snakemake.params.underground_storage[
        "energy_density_MWh_per_m3"
    ]
    cavern_height_bedded_m = snakemake.params.underground_storage[
        "cavern_height_bedded_m"
    ]
    cavern_diameter_bedded_m = snakemake.params.underground_storage[
        "cavern_diameter_bedded_m"
    ]
    cavern_height_dome_m = snakemake.params.underground_storage["cavern_height_dome_m"]
    cavern_diameter_dome_m = snakemake.params.underground_storage[
        "cavern_diameter_dome_m"
    ]

    rho_min = snakemake.params.underground_storage["rho_min"]
    rho_max = snakemake.params.underground_storage["rho_max"]
    theta_safety = snakemake.params.underground_storage["theta_safety"]
    lhv_h2 = snakemake.params.underground_storage["lhv_h2"]

    # Calculate GWh per km² for each cavern type
    if method == "surface":
        gwh_per_km2 = {
            "bedded": compute_gwh_per_km2(
                diameter_m=cavern_diameter_bedded_m,
                height_m=cavern_height_bedded_m,
                energy_density_MWh_per_m3=energy_density_MWh_per_m3,
                surface_area_per_cavern_km2=min_area_km2,
            ),
            "dome": compute_gwh_per_km2(
                diameter_m=cavern_diameter_dome_m,
                height_m=cavern_height_dome_m,
                energy_density_MWh_per_m3=energy_density_MWh_per_m3,
                surface_area_per_cavern_km2=min_area_km2,
            ),
        }
        filtered["gwh_per_km2"] = filtered["salt_type"].map(gwh_per_km2)
        filtered["capacity_gwh"] = filtered["Area_km2"] * filtered["gwh_per_km2"]

    if method == "physics":
        capacity_per_cavern = {
            "bedded": compute_physical_capacity(
                cavern_diameter_bedded_m,
                cavern_height_bedded_m,
                rho_min,
                rho_max,
                theta_safety,
                lhv_h2,
            ),
            "dome": compute_physical_capacity(
                cavern_diameter_dome_m,
                cavern_height_dome_m,
                rho_min,
                rho_max,
                theta_safety,
                lhv_h2,
            ),
        }
        separation_distance_m = {
            "bedded": 4 * cavern_diameter_bedded_m,
            "dome": 4 * cavern_diameter_dome_m,
        }
        caverns_per_km2 = {
            salt_type: 1e6 / (sep**2)  # km² to m², then divide by area per cavern
            for salt_type, sep in separation_distance_m.items()
        }
        # Apply per row
        filtered["capacity_gwh"] = filtered.apply(
            lambda row: row.Area_km2
            * caverns_per_km2[row.salt_type]
            * capacity_per_cavern[row.salt_type],
            axis=1,
        )

    filtered["storage_type"] = "salt_cavern"
    return filtered


def concat_gdf(gdf_list):
    """
    Concatenate multiple geopandas dataframes with common coordinate reference
    system (crs).
    """
    return gpd.GeoDataFrame(pd.concat(gdf_list), crs=geo_crs)


def load_bus_regions(onshore_path, offshore_path):
    """
    Load on- and offshore regions and concat with region_type.
    """
    bus_regions_offshore = gpd.read_file(offshore_path)
    bus_regions_offshore["region_type"] = "offshore"

    bus_regions_onshore = gpd.read_file(onshore_path)
    bus_regions_onshore["region_type"] = "onshore"

    bus_regions = concat_gdf([bus_regions_offshore, bus_regions_onshore])
    bus_regions = bus_regions.dissolve(by=["name", "region_type"], aggfunc="sum")
    return bus_regions


def area(gdf):
    """
    Returns area of GeoDataFrame geometries in square kilometers.
    """
    return gdf.to_crs(area_crs).area / 1e6  # in km²


def salt_cavern_potential_by_region(cavern, regions):
    # calculate area of caverns shapes
    cavern["area_caverns"] = area(cavern)

    overlay = gpd.overlay(regions.reset_index(), cavern, keep_geom_type=True)

    # calculate share of cavern area inside region
    overlay["share"] = area(overlay) / overlay["area_caverns"]

    overlay["e_nom"] = overlay.eval("capacity_gwh * share * area_caverns / 1000")
    cavern_regions = overlay.pivot_table(
        index="name", columns="region_type", values="e_nom", aggfunc="sum"
    ).fillna(0.0)
    return cavern_regions


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_salt_cavern_potentials", clusters="10", simpl=""
        )

    area_crs = snakemake.params.crs["area_crs"]
    geo_crs = snakemake.params.crs["geo_crs"]
    copernicus_crs = COPERNICUS_CRS
    distance_crs = snakemake.params.crs["distance_crs"]
    paths = snakemake.input

    # Load potash deposits shapefile
    gdf = download_potash_data()

    min_area_km2 = snakemake.params.underground_storage["min_area_km2"]

    fn_onshore = snakemake.input.regions_onshore
    fn_offshore = snakemake.input.regions_offshore

    regions = load_bus_regions(fn_onshore, fn_offshore)

    cavern_potential = estimate_h2_potential_from_potash(
        gdf,
        regions,
        min_area_km2,
    )

    # Compute potential
    cavern_regions = salt_cavern_potential_by_region(cavern_potential, regions)

    to_csv_nafix(cavern_regions, snakemake.output.h2_cavern)
