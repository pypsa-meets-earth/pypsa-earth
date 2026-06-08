# SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
# Build salt cavern potentials for hydrogen storage.

https://dx.doi.org/10.2139/ssrn.6307406

This module computes the technical hydrogen storage potential in underground
salt caverns based on global potash deposit data.

The workflow performs:

- classification of salt deposits
- land-use exclusion filtering
- cavern capacity estimation
- aggregation to PyPSA bus regions

The resulting dataset provides hydrogen storage potentials in **GWh per region**
for use in PyPSA-Earth energy system models.
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
from _helpers import mock_snakemake, to_csv_nafix

def capsule_volume(diameter_m: float, height_m: float) -> float:
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
    V = π (d / 2)² (h − d) + (4/3) π (d / 2)³
    where:
    - d = diameter
    - h = total height
    - (h − d) = cylindrical height
    - two hemispheres together form a full sphere of diameter d
    """
    radius = diameter_m / 2
    h_cyl = height_m - diameter_m
    volume_cylinder = math.pi * radius**2 * h_cyl
    volume_spheres = (4 / 3) * math.pi * radius**3
    return volume_cylinder + volume_spheres


def compute_physical_capacity(
    diameter_m: float,
    height_m: float,
    rho_min: float,
    rho_max: float,
    theta_safety: float,
    lhv_h2: float,
) -> float:
    """
    Compute the physical hydrogen storage energy capacity of an underground cavern.

    The function estimates the usable hydrogen energy stored in a salt cavern
    based on cavern geometry, hydrogen density limits, and the lower heating value.

    Parameters
    ----------
    diameter_m : float
        Diameter of the cavern in meters.
    height_m : float
        Height of the cavern in meters.
    rho_min : float
        Minimum hydrogen density in the cavern (kg/m³).
    rho_max : float
        Maximum hydrogen density in the cavern (kg/m³).
    theta_safety : float
        Usable fraction of the working gas capacity.
    lhv_h2 : float
        Lower heating value (LHV) of hydrogen in MJ/kg.

    Returns
    -------
    float
        Physical energy capacity of the cavern in GWh.

    Notes
    -----
    The cavern geometry is approximated as a capsule-shaped structure.
    The working gas mass is calculated from the difference between maximum
    and minimum hydrogen density multiplied by the cavern volume.
    """
    v_cavern = capsule_volume(diameter_m, height_m)
    m_working = (rho_max - rho_min) * v_cavern * theta_safety
    return m_working * lhv_h2 / 1e6  # GWh


def compute_gwh_per_km2(
    diameter_m: float,
    height_m: float,
    energy_density_MWh_per_m3: float = 0.045,
    surface_area_per_cavern_km2: float = 13.0,
) -> float:
    """
    Compute the hydrogen storage potential in GWh per km² for a capsule-shaped salt cavern.

    Parameters
    ----------
    diameter_m : float
        Diameter of the cavern in meters.
    height_m : float
        Total height of the cavern in meters.
    energy_density_MWh_per_m3 : float, optional
        Usable hydrogen energy per cubic meter (default: 0.045 MWh/m³).
    surface_area_per_cavern_km2 : float, optional
        Surface area required per cavern in km² (default: 13 km²).

    Returns
    -------
    float
        Hydrogen storage potential in GWh per km² of surface area.
    """
    volume_m3 = capsule_volume(diameter_m, height_m)
    energy_per_cavern_MWh = volume_m3 * energy_density_MWh_per_m3
    caverns_per_km2 = 1 / surface_area_per_cavern_km2
    return (energy_per_cavern_MWh * caverns_per_km2) / 1000  # Convert to GWh


def classify_salt_type(gdf: gpd.GeoDataFrame) -> pd.Series:
    """
    Classify the type of salt deposit based on geological and deposit information.

    The function adds a new column `salt_type` to the GeoDataFrame with one of
    the following categories:
    - excluded_brine: if the `Dep_type` column contains the word "brine"
    - dome: if the `Dep_type` or `Geology` column contains the word "halokinetic"
    - bedded: if the `Dep_type` or `Geology` column contains the words
      "stratabound" or "evaporite"
    - unknown: if none of the above conditions apply

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        GeoPandas dataframe containing at least the columns `Dep_type` and `Geology`.

    Returns
    -------
    pd.Series
        Series containing the classified salt deposit types.
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


def apply_landuse_exclusions(
    gdf: gpd.GeoDataFrame, regions: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Apply land-use exclusion zones to underground salt storage areas.

    Regions intersecting with specified Copernicus land cover classes
    are removed from the candidate storage areas.

    Workflow:
    1. Apply a negative buffer to each salt storage geometry depending on salt type.
    2. Load Copernicus land cover raster data.
    3. Extract raster cells corresponding to exclusion land-use codes.
    4. Convert raster cells to vector geometries.
    5. Merge exclusion geometries and subtract them from the salt areas.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        Salt storage candidate areas containing a `salt_type` column.
    regions : gpd.GeoDataFrame
        Regions used to clip the Copernicus raster.

    Returns
    -------
    gpd.GeoDataFrame
        Filtered GeoDataFrame excluding areas intersecting with Copernicus
        land-use exclusion zones.

    Notes
    -----
    Buffer distances and exclusion grid codes are defined in
    `snakemake.params.underground_storage`.
    """
    cop = snakemake.params.underground_storage["copernicus"]
    salt_buffer_m_bedded = snakemake.params.underground_storage["salt_buffer_m_bedded"]
    salt_buffer_m_dome = snakemake.params.underground_storage["salt_buffer_m_dome"]

    buffer_map = cop.get("distance_buffer_by_code", {}) or {}
    valid_codes = set(cop["grid_codes"])

    gdf = gdf.to_crs(distance_crs).copy()

    # Apply buffer by salt type before masking
    gdf["geometry"] = gdf.apply(
        lambda row: row.geometry.buffer(
            -salt_buffer_m_bedded if row.salt_type == "bedded" else -salt_buffer_m_dome
        ),
        axis=1,
    )
    gdf = gdf[gdf.is_valid & ~gdf.is_empty]

    # Loading Copernicus raster
    with rioxarray.open_rasterio(paths.copernicus, masked=True) as rds:
        rds_clip = rds.squeeze().rio.clip_box(*regions.total_bounds)
        rds_reproj = rds_clip.rio.reproject(distance_crs)
        clc = rds_reproj.load()
        transform = clc.rio.transform()

        # Combined mask for all relevant grid codes
        mask = np.isin(clc.data, list(valid_codes))

        if not np.any(mask):
            return gdf.to_crs(geo_crs)

        # Single vectorization step, preserving raster values
        shapes_gen = rasterio.features.shapes(clc.data, mask=mask, transform=transform)

        exclusion_geoms = []

        for geom, value in shapes_gen:
            if value not in valid_codes:
                continue

            g = shapely.geometry.shape(geom)

            # Apply code-specific buffer
            buffer_dist = buffer_map.get(str(value), 0)
            if buffer_dist > 0:
                g = g.buffer(buffer_dist)

            if g.is_valid and not g.is_empty:
                exclusion_geoms.append(g)

    if not exclusion_geoms:
        return gdf.to_crs(geo_crs)

    # Merge all exclusion geometries
    union_geom = shapely.ops.unary_union(exclusion_geoms)

    # Subtract exclusions from salt geometries
    gdf["geometry"] = gdf.geometry.difference(union_geom)
    gdf = gdf[gdf.is_valid & ~gdf.is_empty]

    return gdf.to_crs(geo_crs)


def estimate_h2_potential_from_potash(
    potash_gdf: gpd.GeoDataFrame, regions: gpd.GeoDataFrame, min_area_km2: float = 13.0
) -> gpd.GeoDataFrame:
    """
    Estimate technical hydrogen storage potential from potash deposits.

    The function filters and processes potash tract data to determine
    potential salt cavern hydrogen storage capacity.

    Processing Steps:
    - Classify salt deposit types.
    - Filter valid salt formations (bedded or dome).
    - Apply land-use exclusions.
    - Remove small remaining structures.
    - Estimate storage potential using either a surface-based or physics-based method.

    Parameters
    ----------
    potash_gdf : gpd.GeoDataFrame
        Raw potash GIS dataset.
    regions : gpd.GeoDataFrame
        Bus regions used for spatial aggregation.
    min_area_km2 : float
        Minimum salt structure area required to be considered.

    Returns
    -------
    gpd.GeoDataFrame
        Filtered and annotated GeoDataFrame including estimated hydrogen
        storage capacity (`capacity_gwh`).
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


def concat_gdf(gdf_list: list) -> gpd.GeoDataFrame:
    """
    Concatenate multiple GeoDataFrames with a shared CRS.

    Parameters
    ----------
    gdf_list : list
        List of GeoDataFrames.

    Returns
    -------
    gpd.GeoDataFrame
        Combined GeoDataFrame with the global `geo_crs`.
    """
    return gpd.GeoDataFrame(pd.concat(gdf_list), crs=geo_crs)


def load_bus_regions(onshore_path: str, offshore_path: str) -> gpd.GeoDataFrame:
    """
    Load and merge onshore and offshore PyPSA bus regions.

    Parameters
    ----------
    onshore_path : str
        Path to the onshore regions shapefile.
    offshore_path : str
        Path to the offshore regions shapefile.

    Returns
    -------
    gpd.GeoDataFrame
        Combined GeoDataFrame containing both region types.
    """
    bus_regions_offshore = gpd.read_file(offshore_path)
    bus_regions_offshore["region_type"] = "offshore"

    bus_regions_onshore = gpd.read_file(onshore_path)
    bus_regions_onshore["region_type"] = "onshore"

    bus_regions = concat_gdf([bus_regions_offshore, bus_regions_onshore])
    bus_regions = bus_regions.dissolve(by=["name", "region_type"], aggfunc="sum")
    return bus_regions


def area(gdf: gpd.GeoDataFrame) -> pd.Series:
    """
    Compute the area of geometries in a GeoDataFrame.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        Input GeoDataFrame.

    Returns
    -------
    pd.Series
        Area values in square kilometers.
    """
    return gdf.to_crs(area_crs).area / 1e6  # in km²


def salt_cavern_potential_by_region(
    cavern: gpd.GeoDataFrame, regions: gpd.GeoDataFrame
) -> pd.DataFrame:
    """
    Aggregate salt cavern hydrogen storage potentials by region.

    The function overlays cavern storage areas with PyPSA bus regions
    and distributes the storage potential proportionally based on the
    spatial overlap.

    Parameters
    ----------
    cavern : gpd.GeoDataFrame
        Salt cavern storage areas including capacity estimates.
    regions : gpd.GeoDataFrame
        PyPSA bus regions.

    Returns
    -------
    pd.DataFrame
        Regional hydrogen storage potential aggregated by onshore and
        offshore regions.
    """
    # calculate area of caverns shapes
    cavern["area_caverns"] = area(cavern)

    eta_tech = snakemake.params.underground_storage.get("technical_realization_factor", 3.5e-5)

    overlay = gpd.overlay(regions.reset_index(), cavern, keep_geom_type=True)

    # calculate share of cavern area inside region
    overlay["share"] = area(overlay) / overlay["area_caverns"]

    overlay["e_nom"] = overlay.eval("capacity_gwh * share * area_caverns * eta_tech / 1000")
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
    paths = snakemake.input

    # Load potash deposits shapefile
    gdf = gpd.read_file(snakemake.input.potash_shp)

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
