# -*- coding: utf-8 -*-
"""
Build salt cavern potentials for hydrogen storage.

https://doi.org/10.3133/sir20105090S
"""

import os
import shutil
import zipfile
from pathlib import Path

import geopandas as gpd
import pandas as pd
import requests
from _helpers import mock_snakemake


def download_potash_data():
    # URL of the Potash GIS data
    url = "https://pubs.usgs.gov/sir/2010/5090/s/PotashGIS.zip"

    # Download directory
    download_dir = "resources/potash_gis"
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

    # Load the shapefile
    gdf = gpd.read_file(shp_path)

    # Optionally: remove all files except the ones in the same folder as the shapefile
    for item in os.listdir(download_dir):
        full_path = os.path.join(download_dir, item)
        if full_path != os.path.dirname(shp_path):
            if os.path.isdir(full_path):
                shutil.rmtree(full_path)
            else:
                os.remove(full_path)

    return gdf


def estimate_h2_potential_from_potash(
    potash_gdf,
    gwh_per_km2=2.0,
    min_area_km2=5.0,
):
    """
    Estimate underground hydrogen storage potential from Potash tracts.
    """

    halite_keywords = ["halite", "rock salt", "nacl", "salt", "evaporite"]

    # Textbasierte Filterung nach Halit/Salt-Merkmalen in relevanten Spalten
    def is_salt_related(row):
        text = " ".join(
            str(row[col]).lower()
            for col in ["Commodity", "Dep_type", "K_minerals", "Geology"]
            if pd.notnull(row[col])
        )
        return any(kw in text for kw in halite_keywords)

    # Filter nach Geologie und Mindestfläche
    filtered = potash_gdf[
        potash_gdf.apply(is_salt_related, axis=1)
        & (potash_gdf["Area_km2"] > min_area_km2)
    ].copy()

    # Potenzialabschätzung
    filtered["capacity_per_area"] = filtered["Area_km2"] * gwh_per_km2
    filtered["storage_type"] = "salt_cavern"

    return filtered


def concat_gdf(gdf_list, crs="EPSG:4326"):
    """
    Concatenate multiple geopandas dataframes with common coordinate reference
    system (crs).
    """
    return gpd.GeoDataFrame(pd.concat(gdf_list), crs=crs)


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
    return gdf.to_crs(epsg=3035).area / 1e6  # in km²


def salt_cavern_potential_by_region(cavern, regions):
    # calculate area of caverns shapes
    cavern["area_caverns"] = area(cavern)

    overlay = gpd.overlay(regions.reset_index(), cavern, keep_geom_type=True)

    # calculate share of cavern area inside region
    overlay["share"] = area(overlay) / overlay["area_caverns"]

    overlay["e_nom"] = overlay.eval("capacity_per_area * share * area_caverns / 1e6")
    cavern_regions = overlay.pivot_table(
        index="name", columns="region_type", values="e_nom", aggfunc="sum"
    ).fillna(0.0)
    return cavern_regions


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_salt_cavern_potentials", clusters="20", simpl=""
        )

    # Load potash deposits shapefile
    gdf = download_potash_data()

    cavern_potential = estimate_h2_potential_from_potash(
        gdf,
        gwh_per_km2=2.0,
        min_area_km2=5.0,
    )

    fn_onshore = snakemake.input.regions_onshore
    fn_offshore = snakemake.input.regions_offshore

    regions = load_bus_regions(fn_onshore, fn_offshore)

    # Compute potential
    cavern_regions = salt_cavern_potential_by_region(cavern_potential, regions)

    cavern_regions.to_csv(snakemake.output.h2_cavern)
