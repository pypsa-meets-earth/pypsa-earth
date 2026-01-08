# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
# -*- coding: utf-8 -*-
"""
This script handles the downloading and processing of global building data.
"""
import os

import country_converter as coco
import geopandas as gpd
import pandas as pd
from _helpers import (
    BASE_DIR,
    configure_logging,
    create_logger,
)
from shapely import geometry
from tqdm import tqdm

cc = coco.CountryConverter()

logger = create_logger(__name__)


def download_global_buildings_url(update=False):
    """
    Downloads or retrieves the global building URLs from a CSV file, specifically
    from the Microsoft Global Buildings dataset (https://github.com/microsoft/GlobalMLBuildingFootprints).

    Parameters
    ----------
    update : bool, optional
        If True, forces a re-download of the URL list. The default is False.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the global building URLs and country codes.
    """
    global_buildings_path = os.path.join(
        BASE_DIR,
        "data",
        "global_buildings",
        "global_buildings_url.csv",
    )

    if not os.path.exists(global_buildings_path) or update is True:

        country_name = {
            "CongoDRC": "Democratic Republic of the Congo",
            "RepublicoftheCongo": "Republic of the Congo",
            "FYROMakedonija": "North Macedonia",
            "KingdomofSaudiArabia": "Saudi Arabia",
            "SultanateofOman": "Oman",
            "IsleofMan": "Isle of Man",
        }

        logger.info("Downloading Global Buildings URL")
        df_url = pd.read_csv(
            "https://minedbuildings.z5.web.core.windows.net/global-buildings/dataset-links.csv",
            dtype=str,
        )
        df_url["Country"] = [
            cc.convert(names=country, to="ISO2")
            for country in df_url["Location"].replace(country_name)
        ]
        df_url = df_url[df_url["Country"] != "not found"]

        # create global_buildings directory
        os.makedirs(os.path.dirname(global_buildings_path), exist_ok=True)

        df_url.to_csv(global_buildings_path)

    else:
        df_url = pd.read_csv(global_buildings_path, index_col=0)

    return df_url


def get_building_area_center(df, crs):
    """
    Calculates the area and centroid of buildings from a DataFrame of building geometries.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing building URLs and other information.
    crs : dict
        A dictionary containing coordinate reference systems (CRS) for 'geo_crs', 'area_crs', and 'distance_crs'.

    Yields
    ------
    pandas.DataFrame
        A DataFrame with 'area', 'x' and 'y' based on the centroid for each building.
    """
    tqdm_kwargs = dict(
        ascii=False,
        unit=" quadrants",
        desc="Merge Buildings ",
    )

    for i in tqdm(df.index, **tqdm_kwargs):
        geo_df = pd.read_json(df.loc[i, "Url"], lines=True)
        gdf = gpd.GeoDataFrame(
            geometry=geo_df["geometry"].apply(geometry.shape), crs=crs["geo_crs"]
        )
        area = gdf.to_crs(crs["area_crs"]).geometry.area.astype(int).to_list()
        center = (
            gdf.to_crs(crs["distance_crs"])
            .geometry.centroid.to_crs(crs["geo_crs"])
            .to_list()
        )

        yield pd.DataFrame(
            {"area": area, "x": [c.x for c in center], "y": [c.y for c in center]}
        )


def download_global_buildings(country_code, country_buildings_fn, crs, update=False):
    """
    Downloads global building data for a specific country using links from the
    Microsoft Global Buildings dataset (https://github.com/microsoft/GlobalMLBuildingFootprints).
    The shapes are simplified to a rounded area size and center coordinates,
    and then saved into parquet files. This process is done to optimize computer
    memory usage and speed up computational time.

    Parameters
    ----------
    country_code : str
        The ISO2 country code for which to download building data.
    country_buildings_fn : str
        The file path where the processed building data will be saved.
    crs : dict
        A dictionary containing coordinate reference systems (CRS) for 'geo_crs', 'area_crs', and 'distance_crs'.
    update : bool, optional
        If True, forces a re-download of the URL list. The default is False.
    """
    logger.info(f"Downloading Global Buildings for {country_code}")

    df_url = download_global_buildings_url(update=update)
    df_url = df_url[df_url.Country == country_code]
    df = pd.concat(list(get_building_area_center(df_url, crs)), ignore_index=True)

    logger.info(f"Saving Global Buildings for {country_code}")
    df.to_parquet(country_buildings_fn)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("download_global_buildings", simpl="", clusters=10)

    configure_logging(snakemake)

    crs = snakemake.params.crs
    country = snakemake.wildcards.country
    output = snakemake.output[0]

    download_global_buildings(country, output, crs, update=False)
