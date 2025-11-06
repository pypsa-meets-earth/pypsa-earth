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


def download_global_buildings(country_code, geo_crs="EPSG:4326", update=False):
    country_buildings_fn = os.path.join(
        BASE_DIR,
        "data",
        "global_buildings",
        country_code + "_global_buildings_raw.geojson",
    )

    if not os.path.exists(country_buildings_fn) or update is True:
        logger.info(f"Downloading Global Buildings for {country_code}")

        df_url = download_global_buildings_url(update=update)
        df_url = df_url[df_url.Country == country_code]

        tqdm_kwargs = dict(
            ascii=False,
            unit=" quadrants",
            desc="Merge Buildings ",
        )

        geometry_lists = (
            pd.read_json(url, lines=True)["geometry"].apply(geometry.shape).tolist()
            for url in tqdm(df_url["Url"], **tqdm_kwargs)
        )
        geometry_list = list(chain.from_iterable(geometry_lists))

        gdf = gpd.GeoDataFrame(geometry=geometry_list, crs=geo_crs)

        logger.info(f"Saving Global Buildings for {country_code}")
        save_to_geojson(gdf, country_buildings_fn)
    else:
        logger.info(f"Reading Global Buildings for {country_code}")
        gdf = read_geojson(country_buildings_fn)

    return gdf


def calculate_solar_rooftop_area(gb_gdf, shapes, crs, install_ratio, tolerance=100):

    distance_crs = crs["distance_crs"]
    area_crs = crs["area_crs"]

    keys = np.array(sorted(install_ratio.keys()))
    values = np.array([install_ratio[k] for k in keys])

    def get_ratio(area):
        # find the largest key <= area
        idx = np.searchsorted(keys, area, side="right") - 1
        if idx < 0:
            return np.nan  # or 0 if you prefer a default
        return values[idx]

    solar_rooftop_layout = pd.DataFrame()

    for country_code in gb_gdf.keys():
        gdf = gb_gdf[country_code]

        gdf["area"] = gdf.to_crs(area_crs).geometry.area
        gdf["install_ratio"] = gdf["area"].apply(get_ratio)
        gdf["usefull_area"] = gdf["area"] * gdf["install_ratio"]

        gdf.geometry = gdf.to_crs(distance_crs).geometry.centroid

        shapes_country = shapes[shapes.country == country_code].copy()
        shapes_country = shapes_country.to_crs(distance_crs)

        joined = gpd.sjoin(gdf, shapes_country, how="left", predicate="intersects")
        unmatched = joined[joined.name.isna()].copy()
        unmatched = unmatched.drop(["name", "country"], axis=1)

        if not unmatched.empty:
            nearest = gpd.sjoin_nearest(
                unmatched,
                shapes_country,
                max_distance=tolerance * 1e3,
            )

            # Replace the unmatched rows in `joined` with the nearest results
            matched = joined[joined.name.notna()]
            joined = pd.concat([matched, nearest], ignore_index=True)

        usefull_area = joined.groupby("name")["usefull_area"].sum()
        shapes_country["usefull_area"] = shapes_country.index.map(usefull_area)
        solar_rooftop_layout = pd.concat(
            [solar_rooftop_layout, shapes_country["usefull_area"]]
        )

    return solar_rooftop_layout


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("download_global_buildings", simpl="", clusters=10)

    configure_logging(snakemake)
    crs = snakemake.params.crs

    shapes = gpd.read_file(snakemake.input.regions_onshore).set_index("name")[
        ["country", "geometry"]
    ]

    gb_gdf = {
        country_code: download_global_buildings(
            country_code, geo_crs=crs["geo_crs"], update=False
        )
        for country_code in snakemake.params.countries
    }

    solar_rooftop_layout = calculate_solar_rooftop_area(
        gb_gdf,
        shapes,
        crs,
        snakemake.params.install_ratio,
        snakemake.params.tolerance,
    )

    solar_rooftop_layout.to_csv(snakemake.output.solar_rooftop_layout)
