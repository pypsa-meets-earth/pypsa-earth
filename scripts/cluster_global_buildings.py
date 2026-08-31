# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
# -*- coding: utf-8 -*-
"""
Calculate usable solar rooftop area per region from global building data.

This script processes global building footprint data in Parquet format and
aggregates the estimated usable rooftop area for each onshore bus region in a
specified country.

To process large datasets memory-efficiently:
1. The building Parquet file is read sequentially in batches.
2. Useful rooftop areas are calculated via a pre-compiled lookup table.
3. Buildings are mapped to onshore region geometries using spatial joins
   (point-in-polygon with a nearest-neighbor fallback for near-border points).
4. Results are accumulated incrementally per region, keeping memory consumption low.

Inputs
------
snakemake.input.country_buildings : str
    Path to the Parquet file containing building data with columns 'area', 'x', and 'y'.
snakemake.input.regions_onshore : str
    Path to the vector dataset (GeoJSON/GPKG/SHP) containing onshore regions.

Outputs
-------
snakemake.output.solar_rooftop_layout : str
    Path to the output CSV file containing aggregated useful rooftop area per region.

Parameters
----------
snakemake.params.crs : dict
    Dictionary with keys 'distance_crs' (projected CRS) and 'geo_crs' (geographic CRS).
snakemake.params.install_ratio : dict
    Mapping of building area thresholds to usable rooftop fractions.
snakemake.params.solar_rooftop_enable : bool
    Flag indicating whether to execute the calculation or output an empty CSV.
snakemake.params.tolerance : float
    Maximum distance in km to assign unmatched buildings to the nearest region.
GLOBAL_BUILDINGS_BATCH_SIZE : env variable, optional
    Number of rows per Parquet batch (default: 250,000).
"""

import gc
import os
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from _helpers import configure_logging, create_logger

logger = create_logger(__name__)


def _prepare_install_ratio_lookup(install_ratio):
    """
    Prepare sorted threshold arrays for vectorized lookup of install ratios.

    Converts the input dictionary into two synchronized 1D NumPy arrays sorted
    by threshold key, enabling fast 1D array operations (e.g. via ``np.searchsorted``).

    Parameters
    ----------
    install_ratio : dict
        Mapping of capacity/threshold keys to their corresponding install ratios
        (keys and values must be convertible to float).

    Returns
    -------
    keys : numpy.ndarray
        Sorted 1D array of float64 keys/thresholds.
    values : numpy.ndarray
        1D array of float64 values corresponding to the sorted keys.
    """
    install_ratio = {float(k): float(v) for k, v in install_ratio.items()}
    keys = np.array(sorted(install_ratio.keys()), dtype="float64")
    values = np.array([install_ratio[k] for k in keys], dtype="float64")
    return keys, values


def _calculate_useful_area(area, keys, values):
    """
    Calculate useful rooftop area from building footprint area.

    The largest threshold <= building area is used, matching the logic of the
    original get_ratio() implementation. Vectorized using binary search on the
    pre-sorted lookup arrays.

    Parameters
    ----------
    area : array-like
        Building footprint area values (convertible to a float64 NumPy array).
    keys : numpy.ndarray
        Sorted 1D array of threshold keys, as prepared by
        ``_prepare_install_ratio_lookup``.
    values : numpy.ndarray
        1D array of install ratios corresponding to ``keys``.

    Returns
    -------
    numpy.ndarray
        Array of calculated useful rooftop areas (in the same unit as input area).
    """
    area = np.asarray(area, dtype="float64")
    idx = np.searchsorted(keys, area, side="right") - 1

    ratios = np.zeros(len(area), dtype="float64")
    valid = idx >= 0
    ratios[valid] = values[idx[valid]]

    return area * ratios


def _prepare_shapes_country(shapes, country_code, distance_crs):
    """
    Prepare country-specific region geometries for spatial joins.

    Filters the input region shapes by country code, ensures the required
    'name' column exists, and reprojects the geometries to the target coordinate
    reference system (CRS).

    Required columns
    ----------------
    country : str
        Country code used to filter the regions.
    geometry : geometry
        Geospatial geometries of the regions.
    name : str
        Region identifiers (can also be present as the index).

    Parameters
    ----------
    shapes : geopandas.GeoDataFrame
        GeoDataFrame containing region shapes for one or multiple countries.
    country_code : str
        Country code to filter for (e.g. 'DE', 'FR').
    distance_crs : str, int, or pyproj.CRS
        Coordinate reference system to project the output geometries into.

    Returns
    -------
    geopandas.GeoDataFrame
        Filtered and reprojected GeoDataFrame containing only the 'name',
        'country', and 'geometry' columns for the specified country.

    Raises
    ------
    ValueError
        If no regions are found for the given country code.
    KeyError
        If no 'name' column or index exists in the input GeoDataFrame.
    """
    shapes_country = shapes[shapes.country == country_code].copy()

    if shapes_country.empty:
        raise ValueError(
            f"No onshore bus regions found for country '{country_code}'. "
            "Check the 'country' column in regions_onshore."
        )

    # The original script uses set_index("name"). For sjoin we need 'name'
    # as a normal column, and for output we need it as the series index.
    if "name" not in shapes_country.columns:
        shapes_country = shapes_country.reset_index()

    if "name" not in shapes_country.columns:
        raise KeyError("The regions_onshore file must contain a 'name' column or index.")

    shapes_country = shapes_country[["name", "country", "geometry"]].to_crs(distance_crs)
    return shapes_country


def _add_grouped_area(accumulator, joined):
 """
    Add grouped useful_area values from one spatial join result to accumulator.

    Groups the spatial join results by region name, sums the useful rooftop area
    per region, and adds the result element-wise to the existing accumulator Series.

    Required columns
    ----------------
    name : str
        Region identifier used for grouping.
    usefull_area : float
        Calculated useful rooftop area values to sum up.

    Parameters
    ----------
    accumulator : pandas.Series
        Series mapping region names to their accumulated useful rooftop area.
    joined : geopandas.GeoDataFrame or pandas.DataFrame
        Result of a spatial join containing matched building/grid features and regions.

    Returns
    -------
    pandas.Series
        Updated accumulator Series with added useful rooftop areas.
    """
    if joined.empty or "name" not in joined.columns:
        return accumulator

    joined = joined[joined["name"].notna()]
    if joined.empty:
        return accumulator

    grouped = joined.groupby("name")["usefull_area"].sum()
    return accumulator.add(grouped, fill_value=0.0)


def calculate_solar_rooftop_area(
    country_buildings,
    country_code,
    shapes,
    output,
    crs,
    install_ratio,
    tolerance=100,
    batch_size=250_000,
):
    """
    Calculate usable solar rooftop area per region from global building data.

    This memory-efficient implementation avoids loading all buildings into RAM.
    It reads the parquet file in batches, performs the spatial join batch-wise,
    and keeps only one aggregated Series in memory.
    """
    distance_crs = crs["distance_crs"]
    geo_crs = crs["geo_crs"]

    keys, values = _prepare_install_ratio_lookup(install_ratio)

    shapes_country = _prepare_shapes_country(shapes, country_code, distance_crs)

    accumulator = pd.Series(
        0.0,
        index=pd.Index(shapes_country["name"].to_numpy(), name="name"),
        name="usefull_area",
        dtype="float64",
    )

    parquet_file = pq.ParquetFile(country_buildings)

    required_columns = ["area", "x", "y"]
    available_columns = set(parquet_file.schema.names)
    missing_columns = [c for c in required_columns if c not in available_columns]
    if missing_columns:
        raise KeyError(
            f"Missing columns in {country_buildings}: {missing_columns}. "
            "Expected columns are 'area', 'x', and 'y'."
        )

    logger.info(
        f"Aggregating global buildings for {country_code} in batches of "
        f"{batch_size:,} rows"
    )

    for batch_no, batch in enumerate(
        parquet_file.iter_batches(columns=required_columns, batch_size=batch_size),
        start=1,
    ):
        df = batch.to_pandas()

        if df.empty:
            continue

        df = df.dropna(subset=["area", "x", "y"])
        if df.empty:
            continue

        df["usefull_area"] = _calculate_useful_area(df["area"].to_numpy(), keys, values)

        # Buildings with zero useful area cannot contribute to the result.
        df = df[df["usefull_area"] > 0.0]
        if df.empty:
            continue

        # Create point geometries only for the current batch.
        df = df.reset_index(drop=True)
        df["building_id"] = np.arange(len(df))

        gdf = gpd.GeoDataFrame(
            df[["building_id", "usefull_area"]],
            geometry=gpd.points_from_xy(df["x"], df["y"]),
            crs=geo_crs,
        ).to_crs(distance_crs)

        # First try exact point-in-polygon assignment.
        joined = gpd.sjoin(
            gdf,
            shapes_country[["name", "country", "geometry"]],
            how="left",
            predicate="intersects",
        )
        joined = joined.sort_values("building_id").drop_duplicates("building_id")
        matched = joined[joined["name"].notna()]
        accumulator = _add_grouped_area(accumulator, matched)

        # For buildings not inside a region, use nearest region within tolerance.
        unmatched = joined[joined["name"].isna()]
        if not unmatched.empty:
            cols_to_drop = [
                c for c in ["name", "country", "index_right"] if c in unmatched.columns
            ]
            unmatched = unmatched.drop(columns=cols_to_drop)

            nearest = gpd.sjoin_nearest(
                unmatched,
                shapes_country[["name", "country", "geometry"]],
                how="left",
                max_distance=tolerance * 1e3,
            )
            nearest = nearest.sort_values("building_id").drop_duplicates("building_id")

            accumulator = _add_grouped_area(accumulator, nearest)

            del nearest, unmatched

        if batch_no % 10 == 0:
            logger.info(f"Processed {batch_no} parquet batches")

        del df, gdf, joined, matched, batch
        gc.collect()

    Path(output).parent.mkdir(parents=True, exist_ok=True)

    # Same output structure as the original:
    # CSV with index name 'name' and column 'usefull_area'.
    accumulator = accumulator.reindex(shapes_country["name"].to_numpy()).fillna(0.0)
    accumulator.index.name = "name"
    accumulator.name = "usefull_area"
    accumulator.to_csv(output)

    logger.info(f"Saved solar rooftop layout to {output}")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "cluster_global_buildings",
            simpl="",
            clusters=10,
            country="NG",
        )

    configure_logging(snakemake)

    country_code = snakemake.wildcards.country

    logger.info(f"Reading Global Buildings for {country_code}")

    shapes = gpd.read_file(snakemake.input.regions_onshore).set_index("name")[
        ["country", "geometry"]
    ]

    if snakemake.params.solar_rooftop_enable:
        logger.info(f"Calculate solar rooftop area for {country_code}")

        # Can be changed without editing the Snakefile, e.g.:
        # GLOBAL_BUILDINGS_BATCH_SIZE=100000 snakemake -j1 ...
        batch_size = int(os.environ.get("GLOBAL_BUILDINGS_BATCH_SIZE", "250000"))

        calculate_solar_rooftop_area(
            snakemake.input.country_buildings,
            country_code,
            shapes,
            snakemake.output.solar_rooftop_layout,
            crs=snakemake.params.crs,
            install_ratio=snakemake.params.install_ratio,
            tolerance=snakemake.params.tolerance,
            batch_size=batch_size,
        )
    else:
        logger.warning("solar_rooftop_enable is False. Writing an empty output file.")
        output = Path(snakemake.output.solar_rooftop_layout)
        output.parent.mkdir(parents=True, exist_ok=True)
        pd.Series(dtype="float64", name="usefull_area").to_csv(output)