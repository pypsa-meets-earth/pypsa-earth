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
)

cc = coco.CountryConverter()

logger = create_logger(__name__)


def calculate_solar_rooftop_area(
    gdf,
    country_code,
    shapes,
    output,
    crs,
    install_ratio,
    tolerance=100,
):
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
    solar_rooftop_layout = shapes_country["usefull_area"]

    # Save file
    solar_rooftop_layout.to_csv(output)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "cluster_global_buildings", simpl="", clusters=10, country="NG"
        )

    configure_logging(snakemake)
    country_code = snakemake.wildcards.country

    # Retrieve files
    logger.info(f"Reading Global Buildings for {country_code}")
    gdf = read_geojson(snakemake.input.country_buildings)
    shapes = gpd.read_file(snakemake.input.regions_onshore).set_index("name")[
        ["country", "geometry"]
    ]

    if snakemake.params.solar_rooftop_enable:
        logger.info(f"Calculate solar rooftop area for {country_code}")

        calculate_solar_rooftop_area(
            gdf,
            country_code,
            shapes,
            snakemake.output.solar_rooftop_layout,
            crs=snakemake.params.crs,
            install_ratio=snakemake.params.install_ratio,
            tolerance=snakemake.params.tolerance,
        )
