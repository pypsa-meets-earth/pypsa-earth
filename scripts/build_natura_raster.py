# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Converts vectordata or known as shapefiles (i.e. used for geopandas/shapely) to our cutout rasters. The `Protected Planet Data <https://www.protectedplanet.net/en/thematic-areas/wdpa?tab=WDPA>`_ on protected areas is aggregated to all cutout regions.

Relevant Settings
-----------------

.. code:: yaml

    renewable:
        {technology}:
            cutout:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`renewable_cf`

Inputs
------

- ``data/landcover/world_protected_areas/*.shp``: shapefiles representing the world protected areas, such as the `World Database of Protected Areas (WDPA) <https://www.protectedplanet.net/en/thematic-areas/wdpa?tab=WDPA>`_.

    .. image:: /img/natura.png
        :width: 33 %

Outputs
-------

- ``resources/natura/natura.tiff``: Rasterized version of the world protected areas, such as `WDPA <https://www.protectedplanet.net/en/thematic-areas/wdpa?tab=WDPA>`_ natural protection areas to reduce computation times.

    .. image:: /img/natura.png
        :width: 33 %

Description
-----------
To operate the script you need all input files.

This script collects all shapefiles available in the folder `data/landcover/*` describing regions of protected areas,
merges them to one shapefile, and create a rasterized version of the region, that covers the region described by the cutout.
The output is a raster file with the name `natura.tiff` in the folder `resources/natura/`.
"""
import os

import atlite
import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio as rio
from _helpers import configure_logging, create_logger
from rasterio.features import geometry_mask, rasterize
from rasterio.warp import transform_bounds
from rasterio.windows import from_bounds
from shapely.ops import unary_union

logger = create_logger(__name__)


CUTOUT_CRS = "EPSG:4326"


def get_relevant_regions(country_shapes, offshore_shapes, natura_crs):
    """
    Merge the country_shapes and the offshore_shapes into one GeoDataFrame.
    Additionally add a buffer to ensure all relevant regions are included.

    Returns
    -------
    regions : GeoDataFrame with a unified "multipolygon"
    """

    # unify the country_shapes and offshore_shapes to select the regions of interest
    # load country shapes
    countries_gdf = gpd.read_file(country_shapes).to_crs(natura_crs)
    countries = countries_gdf.geometry.union_all()

    # load offshore shapes
    offshore_gdf = gpd.read_file(offshore_shapes).to_crs(natura_crs)
    offshore = offshore_gdf.geometry.union_all()

    # combine countries and offshore regions into one merged geometry
    buffer = 1e5  # 100 km, adjust should important regions be missed
    merged_regions = unary_union([countries.buffer(buffer), offshore.buffer(buffer)])

    regions = gpd.GeoDataFrame(geometry=[merged_regions], crs=natura_crs)

    return regions


def get_fileshapes(list_paths, accepted_formats=(".shp",)):
    "Function to parse the list of paths to include shapes included in folders, if any"

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


def determine_region_xXyY(cutout_name, regions, natura_size, out_logging):
    '''
    Determine the bounds of the analyzed regions depending on the natura_size parameter.
    "global" includes the entire world, "cutout" the extend of the cutout, and
    "countries" only includes the bounds of the requested countries and their offshore regions.

    Returns
    -------
    cutout_xXyY : List including the bounds
    '''

    if out_logging:
        logger.info("Stage 1/5: Determine cutout boundaries")

    cutout = atlite.Cutout(cutout_name)
    assert cutout.crs == CUTOUT_CRS
    dx, dy = cutout.dx, cutout.dy

    match natura_size:
        case "global":
            return [-180, 180, -90, 90]
        case "cutout":
            x, X, y, Y = cutout.extent
        case "countries":
            x, y, X, Y = regions.to_crs(CUTOUT_CRS).total_bounds
        case _:
            raise ValueError(f"Provided an unknown value for the parameter 'natura_size': {natura_size}")

    cutout_xXyY = [
        np.clip(x - dx / 2.0, -180, 180),
        np.clip(X + dx / 2.0, -180, 180),
        np.clip(y - dy / 2.0, -90, 90),
        np.clip(Y + dy / 2.0, -90, 90),
    ]
    return cutout_xXyY


def get_transform_and_shape(bounds, res, out_logging):
    if out_logging:
        logger.info("Stage 2/5: Get transform and shape")
    left, bottom = [(b // res) * res for b in bounds[:2]]
    right, top = [(b // res + 1) * res for b in bounds[2:]]
    # "latitude, longitude" coordinates order
    shape = int((top - bottom) // res), int((right - left) // res)
    transform = rio.Affine(res, 0, left, 0, -res, top)
    return transform, shape



def decide_bigtiff_flag(out_shape, dtype="uint8", safety_factor=1.1):
    '''
    Decide whether BIGTIFF should be "YES" or "NO" based on raster shape.
    BIGTIFF is required for filesizes larger than 4 GB.

    Returns
    -------
    str: "YES" if the estimated size is larger than 4 GB, else "NO".
    '''

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

    regions = get_relevant_regions(country_shapes, offshore_shapes, natura_crs)

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

        for file in shp_files:
            shp = gpd.read_file(file).to_crs(natura_crs)
            total_files += 1
            try:
                shp.geometry = shp.geometry.make_valid()
                read_files += 1
                logger.info(f"Successfully read file {file}")
            except Exception as e:
                logger.warning(f"Error reading file {file}: {e}")
                continue

            max_row = out_shape[0] - 1  # height in pixels
            max_col = out_shape[1] - 1  # width in pixels

            for row0 in range(0, max_row, window_size):  # rows = y direction
                row1 = min(row0 + window_size, max_row)
                ymin = top - row1 * res
                ymax = top - row0 * res

                for col0 in range(0, max_col, window_size):  # cols = x direction
                    col1 = min(col0 + window_size, max_col)
                    xmin = left + col0 * res
                    xmax = left + col1 * res

                    window = from_bounds(xmin, ymin, xmax, ymax, transform=transform)

                    window_shp = shp.cx[xmin:xmax, ymin:ymax]

                    if window_shp.empty:
                        continue

                    geom = unary_union(window_shp.geometry)

                    data = dst.read(1, window=window)

                    mask = rasterize(
                        [(geom, 1)],
                        out_shape=data.shape,
                        transform=dst.window_transform(window),
                        fill=0,
                        dtype="uint8"
                    )

                    data[mask == 1] = 1

                    dst.write(data, 1, window=window)

                logger.info(f"Done writing rows {row0}-{row1} out of {max_row}")

            logger.info(f"Done writing file {file}")