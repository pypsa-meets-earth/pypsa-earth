# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
# -*- coding: utf-8 -*-
"""
TEST
"""
import os
from itertools import chain

import country_converter as coco
import geopandas as gpd
import numpy as np
import pandas as pd
from _helpers import (
    BASE_DIR,
    configure_logging,
    create_logger,
    read_geojson,
    save_to_geojson,
)
from shapely import geometry
from tqdm import tqdm

cc = coco.CountryConverter()

logger = create_logger(__name__)


def download_global_buildings_url(update=False):
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

        yield pd.DataFrame({"area": area, "center": center})


def download_global_buildings(country_code, country_buildings_fn, crs, update=False):
    logger.info(f"Downloading Global Buildings for {country_code}")

    df_url = download_global_buildings_url(update=update)
    df_url = df_url[df_url.Country == country_code]
    df = pd.concat(list(get_building_area_center(df_url, crs)), ignore_index=True)

    logger.info(f"Saving Global Buildings for {country_code}")
    df.to_csv(country_buildings_fn)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("download_global_buildings", simpl="", clusters=10)

    configure_logging(snakemake)

    crs = snakemake.params.crs
    country = snakemake.wildcards.country
    output = snakemake.output[0]

    download_global_buildings(country, output, crs, update=False)
