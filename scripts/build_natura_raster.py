# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Builds a rasterized Natura and environmental exclusion mask from vector geometries.
The geometries are sourced from [Protected Planet Data](https://www.protectedplanet.net/en/thematic-areas/wdpa?tab=WDPA)
and filtered for the required regions.

Relevant Settings
-----------------

    countries:

    enable:
        build_natura_raster:
        progress_bar:

    crs:
        area_crs:

    natura:
        natura_size:
        natura_resolution:
        window_size:
        buffer_size:

    renewable:
        {technology}:
            cutout:

Inputs
------

- ``data/landcover/world_protected_areas/*.shp``: Vectorized shapefiles representing the world protected areas from [Protected Planet Data](https://www.protectedplanet.net/en/thematic-areas/wdpa?tab=WDPA).
- ``cutouts/{CDIR}/{cutout}.nc``: Atlite cutout specified in ``renewable: {technology}: cutout`` for each technology. The cutout(s) are used to determine geographic coverage and resolution alignment.
- ``resources/{RDIR}/shapes/country_shapes.geojson``: Onshore country geometries used to determine the spatial extent of the rasterization.
- ``resources/{RDIR}/shapes/offshore_shapes.geojson``: Offshore regions associated with the selected countries.

Outputs
-------

- ``resources/{RDIR}/natura.tiff``: Rasterized version of the world protected areas.

Description
-----------
The rule ``build_natura_raster`` converts large collections of vector-based
environmental exclusion geometries into a rasterized binary mask containing:

- ``1`` for raster cells intersecting exclusion geometries, and
- ``0`` elsewhere.

First, the rule determines the spatial extent of the rasterization which is configured through
``natura: natura_size`` and can be based on:

- the full global extent,
- the combined extent of all configured renewable technology cutouts
  (``renewable: {technology}: cutout``), or
- the combined extent of all selected country (``countries``) and their offshore regions.

In case of the ``countries`` extent, the country and offshore geometries are reprojected into (``crs: area_crs``).
Additionally, a buffer (``natura: buffer_size``) is applied to the geometries to ensure that all regions of interest
are fully included in the rasterization.

The rasterization is performed as follows:

- an empty GeoTIFF is created using the extent determined above and a grid
  resolution defined by ``natura: natura_resolution``,
- the extent is split into smaller windows with a maximum width and height of ``natura: window_size``,
- each input ``.shp`` file is processed independently and rasterized window by window.

This windowed processing strategy reduces memory usage and allows
the rule to run efficiently on standard compute environments.

Notes
-----
This script only runs in the current PyPSA-Earth workflow if ``enable: build_natura_raster`` is ``true``.
Otherwise, the workflow will use the general ``data/natura/natura.tiff`` file, which is copied using the rule ``copy_defaultnatura_tiff``.

There are two main differences between the two options, the data source and the license:

- The rule ``build_natura_raster`` uses data from [Protected Planet Data](https://www.protectedplanet.net/en/thematic-areas/wdpa?tab=WDPA)
  which is generally more up-to-date. However, the [Protected Planet License](https://www.protectedplanet.net/en/legal)
  is less permissive and might not be applicable for your use case.
- The rule ``copy_defaultnatura_tiff`` uses data from [Harvard Dataverse Protected areas (WDPA)](https://doi.org/10.7910/DVN/XIV9BL),
  which was last updated in 2022. This dataset is licensed under the [CC0 1.0 license](https://creativecommons.org/publicdomain/zero/1.0/).
"""

import os

import atlite
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio as rio
from _helpers import configure_logging, create_logger
from rasterio.features import geometry_mask, rasterize
from rasterio.warp import transform_bounds
from rasterio.windows import from_bounds
from shapely.ops import unary_union
from tqdm import tqdm

logger = create_logger(__name__)


CUTOUT_CRS = "EPSG:4326"


def get_relevant_regions(
    country_shapes: str,
    offshore_shapes: str,
    natura_crs: str,
    buffer: float,
) -> gpd.GeoDataFrame:
    """
    Load and merge country and offshore regions into a unified geometry.

    The resulting geometry is buffered to ensure all relevant nearby regions
    are included.

    Parameters
    ----------
    country_shapes : str
        Path to the vector file containing the country geometries.
    offshore_shapes : str
        Path to the vector file containing the offshore geometries.
    natura_crs : str
        Coordinate reference system used for all geometries.
    buffer : float
        Buffer distance applied to both country and offshore geometries.
        Units are determined by ``natura_crs``.

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame containing a single merged geometry representing the
        buffered union of all country and offshore regions.

    Notes
    -----
    Both input datasets are reprojected to ``natura_crs`` before merging.

    Examples
    --------
    ```python
    regions = get_relevant_regions(
        "resources/shapes/country_shapes.geojson",
        "resources/shapes/offshore_shapes.geojson",
        "ESRI:54009",
        10000,
    )
    len(regions)
    # 1
    regions.crs.to_string()
    # "ESRI:54009"
    ```
    """

    # unify the country_shapes and offshore_shapes to select the regions of interest
    # load country shapes
    countries_gdf = gpd.read_file(country_shapes).to_crs(natura_crs)
    countries = countries_gdf.geometry.union_all()

    # load offshore shapes
    offshore_gdf = gpd.read_file(offshore_shapes).to_crs(natura_crs)
    offshore = offshore_gdf.geometry.union_all()

    # combine countries and offshore regions into one merged geometry
    merged_regions = unary_union([countries.buffer(buffer), offshore.buffer(buffer)])

    regions = gpd.GeoDataFrame(geometry=[merged_regions], crs=natura_crs)

    return regions


def get_fileshapes(
    list_paths: list[str],
    accepted_formats: str | tuple[str, ...] = (".shp",),
) -> list[str]:
    """
    Function to parse the list of paths and identify the ones with one of the accepted file formats.

    Parameters
    ----------
    list_paths : list[str]
        List of paths to check.
    accepted_formats : str | tuple[str, ...], optional
        File format(s) to accepted. If a string is provided, only that extension is accepted.
        If a tuple is provided, any of the extensions are accepted. Default is (".shp",).

    Returns
    -------
    list[str]
        List of all paths which are of one of the accepted formats.

    Examples
    --------
    ```python
    paths = ["shape_file.shp", "shape_index.shi"]
    get_fileshapes(paths)
    # ["shape_file.shp"]
    get_fileshapes(paths, ".shi")
    # ["shape_index.shi"]
    ```
    """
    list_fileshapes = []
    for lf in list_paths:
        if os.path.isdir(lf):  # if it is a folder, then list all shapes files contained
            # loop over all dirs and subdirs
            for path, subdirs, files in os.walk(lf):
                # loop over all files
                for subfile in files:
                    # add the subfile if it is a shape file
                    if subfile.endswith(accepted_formats):
                        list_fileshapes.append(os.path.join(path, subfile))

        elif lf.endswith(accepted_formats):
            list_fileshapes.append(lf)

    return list_fileshapes


def determine_region_xXyY(
    cutout_name: str,
    regions: gpd.GeoDataFrame | None,
    natura_size: str,
    out_logging: bool,
) -> list[float]:
    """
    Determine the bounds of the analyzed regions.

    Parameters
    ----------
    cutout_name : str
        Path of the cutout file.
    regions : gpd.GeoDataFrame | None
        GeoDataFrame containing the analyzed regions.
        Required if natura_size is "countries", otherwise ignored.
    natura_size : str
        Flag to determine which region should be used.
        "global" includes the entire world, "cutout" the extent of the cutout, and
        "countries" only includes the bounds of the requested countries and their offshore regions.
    out_logging : bool
        If True, emits progress information via the module logger.

    Returns
    -------
    list[float]
        Bounding box of the region in the format:
        [min_lon, max_lon, min_lat, max_lat].

    Examples
    --------
    ```python
    cutout_path = "cutouts/cutout-2013-era5.nc"
    regions = None

    # Global extent:
    determine_region_xXyY(cutout_path, regions, "global", False)
    # [-180, 180, -90, 90]

    # Cutout extent (Africa):
    determine_region_xXyY(cutout_path, regions, "cutout", False)
    # [-19.8, 67.8, -37.8, 39.6]

    # Countries extent:
    import geopandas as gpd
    from shapely.geometry import Polygon

    regions = gpd.GeoDataFrame(
        geometry=[Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])],
        crs="EPSG:4326",
    )
    determine_region_xXyY(cutout_path, regions, "countries", False)
    # [0.0, 10.0, 0.0, 10.0]
    ```
    """

    if out_logging:
        logger.info("Stage 1/5: Determine cutout boundaries")

    # identify the cell size of the cutout
    cutout = atlite.Cutout(cutout_name)
    assert cutout.crs == CUTOUT_CRS
    dx, dy = cutout.dx, cutout.dy

    # identify the extent of the analyzed region
    match natura_size:
        case "global":
            return [-180, 180, -90, 90]
        case "cutout":
            x, X, y, Y = cutout.extent
        case "countries":
            if regions is None:
                raise ValueError(
                    "Must provide regions when natura_size is 'countries', instead provided None."
                )
            x, y, X, Y = regions.to_crs(CUTOUT_CRS).total_bounds
        case _:
            raise ValueError(
                f"Provided an unknown value for the parameter 'natura_size': {natura_size}"
            )

    # ensure all cells are fully included and that the bounds remain within the latitude and longitude range
    cutout_xXyY = [
        np.clip(x - dx / 2.0, -180, 180),
        np.clip(X + dx / 2.0, -180, 180),
        np.clip(y - dy / 2.0, -90, 90),
        np.clip(Y + dy / 2.0, -90, 90),
    ]
    return cutout_xXyY


def get_transform_and_shape(
    bounds: list[float],
    res: float,
    out_logging: bool,
) -> tuple[rio.Affine, tuple[int, int]]:
    """
    Compute an affine transform and raster shape from spatial bounds and resolution.

    Parameters
    ----------
    bounds : list[float]
        Bounding box in the format:
        [min_lon, min_lat, max_lon, max_lat].
    res : float
        Spatial resolution of the grid (in coordinate units, e.g. degrees).
    out_logging : bool
        If True, emits progress information via the module logger.

    Returns
    -------
    transform : rio.Affine
        Affine transform mapping raster indices (row, col) to spatial coordinates.
    shape : tuple[int, int]
        Raster shape as (n_lat, n_lon), corresponding to (rows, cols).

    Examples
    --------
    ```python
    import rasterio as rio

    # [min_lon, min_lat, max_lon, max_lat]
    bounds = [0.0, 50.0, 3.0, 52.0]
    res = 1.0

    transform, shape = get_transform_and_shape(bounds, res, False)
    shape
    # (2, 3)
    transform
    # Affine(1.0, 0.0, 0.0,
    # 0.0, -1.0, 52.0)
    ```
    """
    if out_logging:
        logger.info("Stage 2/5: Get transform and shape")

    # snap the bounds to the nearest grid
    left, bottom = [(b // res) * res for b in bounds[:2]]
    right, top = [(b // res + 1) * res for b in bounds[2:]]

    # identify the shape of the raster as "latitude, longitude"
    shape = int((top - bottom) // res), int((right - left) // res)
    # prepare the affine transform to map raster indices (row, col) to spatial coordinates
    transform = rio.Affine(res, 0, left, 0, -res, top)

    return transform, shape


def decide_bigtiff_flag(
    out_shape: tuple[int, int],
    dtype: str = "uint8",
    safety_factor: float = 1.1,
) -> str:
    """
    Decide whether a raster requires BIGTIFF storage based on raster shape.

    BIGTIFF is required for GeoTIFF files larger than approximately 4 GB.

    Parameters
    ----------
    out_shape : tuple[int, int]
        Raster shape as (n_rows, n_cols).
    dtype : str, optional
        Data type of the raster values. Default is "uint8".
    safety_factor : float, optional
        Multiplicative buffer applied to the estimated size to account for
        metadata, compression inefficiencies, and file overhead. Default is 1.1.

    Returns
    -------
    str
        "YES" if the estimated size is larger than the BIGTIFF threshold,
        otherwise "NO".

    Notes
    -----
    The size estimate is computed as:

        n_rows * n_cols * bytes_per_pixel * safety_factor

    The BIGTIFF threshold is fixed at 4,000,000,000 bytes (~4 GB).

    Examples
    --------
    ```python
    decide_bigtiff_flag((100, 100), dtype="uint8")
    # "NO"
    ```
    """

    tiff_size = 4_000_000_000

    bytes_per_pixel = np.dtype(dtype).itemsize
    estimated_size = out_shape[0] * out_shape[1] * bytes_per_pixel * safety_factor

    return "YES" if estimated_size > tiff_size else "NO"


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_natura_raster", cutouts=["cutouts/africa-2013-era5.nc"]
        )
    configure_logging(snakemake)

    # get crs
    natura_crs = snakemake.params.area_crs

    out_logging = True
    inputs = snakemake.input
    cutouts = inputs.cutouts
    shapefiles = get_fileshapes(inputs)
    country_shapes = inputs.country_shapes
    offshore_shapes = inputs.offshore_shapes
    natura_size = snakemake.params.natura["natura_size"]
    natura_resolution = snakemake.params.natura["natura_resolution"]
    window_size = snakemake.params.natura["window_size"]
    buffer_size = snakemake.params.natura["buffer_size"]
    disable_progress = snakemake.params.disable_progress

    if natura_size == "countries":
        regions = get_relevant_regions(
            country_shapes, offshore_shapes, natura_crs, buffer_size
        )
    else:
        regions = None

    xs, Xs, ys, Ys = zip(
        *(
            determine_region_xXyY(cutout, regions, natura_size, out_logging=out_logging)
            for cutout in cutouts
        )
    )

    bounds = transform_bounds(
        CUTOUT_CRS, natura_crs, min(xs), min(ys), max(Xs), max(Ys)
    )

    transform, out_shape = get_transform_and_shape(
        bounds, res=natura_resolution, out_logging=out_logging
    )

    res = transform.a
    left = transform.c
    top = transform.f

    # read only .shp snakemake inputs
    shp_files = [string for string in shapefiles if ".shp" in string]
    assert len(shp_files) != 0, "no input shapefiles given"

    # identify if BIGTIFF is require (filesize > 4 GB)
    bigtiff_flag = decide_bigtiff_flag(out_shape)

    # create a 0 filled file
    with rio.open(
        snakemake.output[0],
        "w",
        driver="GTiff",
        dtype=rio.uint8,
        count=1,
        transform=transform,
        crs=natura_crs,
        compress="lzw",
        width=out_shape[1],
        height=out_shape[0],
        nodata=0,
        BIGTIFF=bigtiff_flag,
    ) as dst:
        pass

    with rio.open(snakemake.output[0], "r+") as dst:

        total_files = 0
        read_files = 0

        # read the shapefiles
        for file in tqdm(
            shp_files,
            desc="Processing shapefiles",
            unit="file",
            position=0,
            disable=disable_progress,
        ):
            shp = gpd.read_file(file).to_crs(natura_crs)
            total_files += 1
            try:
                shp.geometry = shp.geometry.make_valid()
                read_files += 1
                logger.info(f"\nSuccessfully read file {file}")
            except Exception as e:
                logger.warning(f"\nError reading file {file}: {e}")
                continue

            max_row = out_shape[0] - 1  # height in pixels
            max_col = out_shape[1] - 1  # width in pixels

            with tqdm(
                total=(max_row // window_size + 1) * (max_col // window_size + 1),
                desc="Rasterizing windows",
                unit="window",
                position=1,
                disable=disable_progress,
            ) as rasterize_pbar:

                # write the raster with a shifting window
                for row0 in range(0, max_row, window_size):  # rows = y direction
                    row1 = min(row0 + window_size, max_row)
                    ymin = top - row1 * res
                    ymax = top - row0 * res

                    for col0 in range(0, max_col, window_size):  # cols = x direction
                        col1 = min(col0 + window_size, max_col)
                        xmin = left + col0 * res
                        xmax = left + col1 * res

                        window = from_bounds(
                            xmin, ymin, xmax, ymax, transform=transform
                        )

                        window_shp = shp.cx[xmin:xmax, ymin:ymax]

                        if window_shp.empty:
                            rasterize_pbar.update(1)
                            continue

                        geom = unary_union(window_shp.geometry)

                        data = dst.read(1, window=window)

                        mask = rasterize(
                            [(geom, 1)],
                            out_shape=data.shape,
                            transform=dst.window_transform(window),
                            fill=0,
                            dtype="uint8",
                        )

                        data[mask == 1] = 1

                        dst.write(data, 1, window=window)

                        rasterize_pbar.update(1)
