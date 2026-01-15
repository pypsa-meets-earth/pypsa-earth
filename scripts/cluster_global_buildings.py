# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
# -*- coding: utf-8 -*-
"""
This script processes global building data to calculate solar rooftop area.
"""
import country_converter as coco
import geopandas as gpd
import numpy as np
import pandas as pd
from _helpers import (
    BASE_DIR,
    configure_logging,
    create_logger,
)
from shapely import Point

cc = coco.CountryConverter()

logger = create_logger(__name__)


def calculate_solar_rooftop_area(
    df,
    country_code,
    shapes,
    output,
    crs,
    install_ratio,
    tolerance=100,
):
    """
    Calculates the usable solar rooftop area for buildings within a country, considering
    that only a portion of each rooftop can be allocated for PV installation based on its size.
    Smaller buildings are less likely to host PV systems, while larger buildings can dedicate
    a greater share of their rooftop to PVs. This approach is based on the methodology
    described in a paper by Hideaki Obane (https://eneken.ieej.or.jp/data/12710.pdf).

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing building areas and center points.
    country_code : str
        The ISO2 country code.
    shapes : geopandas.GeoDataFrame
        GeoDataFrame of regional shapes with 'country' and 'geometry' columns.
    output : str
        The file path to save the calculated solar rooftop area.
    crs : dict
        A dictionary containing coordinate reference systems (CRS) for 'distance_crs' and 'geo_crs'.
    install_ratio : dict
        A dictionary mapping building area thresholds to solar installation ratios.
    tolerance : int, optional
        Maximum distance in kilometers for spatial join to find nearest shapes.
        The default is 100.
    """
    distance_crs = crs["distance_crs"]
    geo_crs = crs["geo_crs"]

    keys = np.array(sorted(install_ratio.keys()))
    values = np.array([install_ratio[k] for k in keys])

    def get_ratio(area):
        # find the largest key <= area
        idx = np.searchsorted(keys, area, side="right") - 1
        if idx < 0:
            return np.nan  # or 0 if you prefer a default
        return values[idx]

    df["usefull_area"] = df["area"] * df["area"].apply(get_ratio)
    gdf = gpd.GeoDataFrame(df, geometry="center", crs=geo_crs).to_crs(distance_crs)
    shapes_country = shapes[shapes.country == country_code].to_crs(distance_crs)

    joined = gpd.sjoin(gdf, shapes_country, how="left", predicate="intersects")
    unmatched = joined[joined.name.isna()]
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

    # Save file
    shapes_country["usefull_area"].to_csv(output)


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
    df = pd.read_parquet(snakemake.input.country_buildings)
    df["center"] = df.apply(lambda row: Point(row["x"], row["y"]), axis=1)

    shapes = gpd.read_file(snakemake.input.regions_onshore).set_index("name")[
        ["country", "geometry"]
    ]

    if snakemake.params.solar_rooftop_enable:
        logger.info(f"Calculate solar rooftop area for {country_code}")

        calculate_solar_rooftop_area(
            df,
            country_code,
            shapes,
            snakemake.output.solar_rooftop_layout,
            crs=snakemake.params.crs,
            install_ratio=snakemake.params.install_ratio,
            tolerance=snakemake.params.tolerance,
        )
