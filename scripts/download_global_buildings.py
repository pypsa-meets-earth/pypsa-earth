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
import gc
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq
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
    Calculate the footprint area and centroid coordinates of buildings.

    Processes Microsoft Global Buildings quadrant files sequentially (one quadrant
    at a time) to keep memory usage low. Converts raw JSON geometries into spatial
    geometries, reprojects them to derive accurate surface area and centroid coordinates,
    and yields the results chunk-by-chunk.

    Required columns
    ----------------
    Url : str
        URL or file path to the Microsoft Global Buildings JSON/GeoJSON-Lines quadrant file.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing links to building quadrant files in a 'Url' column.
    crs : dict
        Dictionary containing CRS definitions with keys 'geo_crs' (geographic CRS,
        e.g., EPSG:4326), 'area_crs' (equal-area projection for area calculation),
        and 'distance_crs' (projected CRS for centroid calculations).

    Yields
    ------
    pandas.DataFrame
        Dataframe for a single quadrant containing:
        - area (int64): Building footprint area.
        - x (float64): Centroid longitude/x-coordinate in geographic CRS.
        - y (float64): Centroid latitude/y-coordinate in geographic CRS.
    """
    tqdm_kwargs = dict(
        ascii=False,
        unit=" quadrants",
        desc="Merge Buildings ",
    )

    for i in tqdm(df.index, **tqdm_kwargs):
        url = df.loc[i, "Url"]

        # Load only one quadrant
        geo_df = pd.read_json(url, lines=True)

        if geo_df.empty:
            del geo_df
            gc.collect()
            continue

        # Convert JSON geometries to shapely geometries
        geometries = geo_df["geometry"].apply(geometry.shape)

        # Free raw JSON dataframe as early as possible
        del geo_df
        gc.collect()

        gdf = gpd.GeoDataFrame(
            geometry=geometries,
            crs=crs["geo_crs"],
        )

        del geometries
        gc.collect()

        # Calculate area in projected CRS
        area_gdf = gdf.to_crs(crs["area_crs"])
        area = area_gdf.geometry.area.round().astype("int64").to_numpy()

        del area_gdf
        gc.collect()

        # Calculate centroids and transform back to geographic CRS
        center = (
            gdf.to_crs(crs["distance_crs"])
            .geometry.centroid
            .to_crs(crs["geo_crs"])
        )

        result = pd.DataFrame(
            {
                "area": area,
                "x": center.x.to_numpy(dtype="float64"),
                "y": center.y.to_numpy(dtype="float64"),
            }
        )

        del gdf, area, center
        gc.collect()

        yield result

        # Free yielded result after it has been written by the caller
        del result
        gc.collect()


def download_global_buildings(country_code, country_buildings_fn, crs, update=False):
    """
    Download global building data for a country and write results incrementally.

    Fetches Microsoft Global Buildings quadrant URLs for the specified country,
    processes each quadrant using ``get_building_area_center``, and streams the
    processed building records directly into a Parquet file to avoid keeping all
    quadrants in RAM at once.

    Parameters
    ----------
    country_code : str
        Country code used to filter building download URLs (e.g. 'NG', 'DE').
    country_buildings_fn : str or path-like
        Output file path where the final Parquet file will be saved.
    crs : dict
        Dictionary containing CRS definitions with keys 'geo_crs', 'area_crs',
        and 'distance_crs'.
    update : bool, optional
        Whether to force updating/re-downloading the building dataset index URL
        file (default is False).

    Returns
    -------
    None

    Raises
    ------
    Exception
        Re-raises any error encountered during processing after cleaning up
        temporary output files.
    """
    logger.info(f"Downloading Global Buildings for {country_code}")

    df_url = download_global_buildings_url(update=update)
    df_url = df_url[df_url.Country == country_code]

    output_path = Path(country_buildings_fn)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    tmp_output_path = output_path.with_suffix(output_path.suffix + ".tmp")

    if tmp_output_path.exists():
        tmp_output_path.unlink()

    if output_path.exists():
        output_path.unlink()

    writer = None
    schema = None
    total_rows = 0

    try:
        for part_df in get_building_area_center(df_url, crs):
            if part_df.empty:
                continue

            # Ensure stable dtypes across all parquet row groups
            part_df = part_df.astype(
                {
                    "area": "int64",
                    "x": "float64",
                    "y": "float64",
                }
            )

            table = pa.Table.from_pandas(part_df, preserve_index=False)

            if writer is None:
                schema = table.schema
                writer = pq.ParquetWriter(
                    tmp_output_path,
                    schema,
                    compression="snappy",
                )
            else:
                table = table.cast(schema, safe=False)

            writer.write_table(table)
            total_rows += table.num_rows

            del part_df, table
            gc.collect()

        if writer is not None:
            writer.close()
            writer = None
        else:
            # Create an empty parquet file if no buildings were found
            empty_schema = pa.schema(
                [
                    ("area", pa.int64()),
                    ("x", pa.float64()),
                    ("y", pa.float64()),
                ]
            )
            empty_table = pa.Table.from_arrays(
                [
                    pa.array([], type=pa.int64()),
                    pa.array([], type=pa.float64()),
                    pa.array([], type=pa.float64()),
                ],
                schema=empty_schema,
            )
            pq.write_table(empty_table, tmp_output_path, compression="snappy")

        tmp_output_path.replace(output_path)

        logger.info(
            f"Saving Global Buildings for {country_code}: "
            f"{total_rows:,} buildings written to {output_path}"
        )

    except Exception:
        if writer is not None:
            writer.close()

        if tmp_output_path.exists():
            tmp_output_path.unlink()

        raise

def get_building_area_center_old(df, crs):
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


def download_global_buildings_old(country_code, country_buildings_fn, crs, update=False):
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
