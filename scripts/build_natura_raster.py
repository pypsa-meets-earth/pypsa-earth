# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors, 2021 PyPSA-Africa Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later
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

- ``data/raw/protected_areas/WDPA_WDOECM_Aug2021_Public_AF_shp-points.shp``: `WDPA <https://en.wikipedia.org/wiki/Natura_2000>`_ World Database for Protected Areas.

    .. image:: ../img/natura.png
        :scale: 33 %

Outputs
-------

- ``resources/natura.tiff``: Rasterized version of `Natura 2000 <https://en.wikipedia.org/wiki/Natura_2000>`_ natural protection areas to reduce computation times.

    .. image:: ../img/natura.png
        :scale: 33 %

Description
-----------
To operate the script you need all input files. The Snakefile describes what goes in and out. Make sure you didn't skip one of these.
Maybe not so obvious is the cutout input. An example is this `africa-2013-era5.nc`

Steps to retrieve the protected area data (as apparently no API is given for the WDPA data):
    - 1. Download the WPDA Dataset: World Database on Protected Areas. UNEP-WCMC and IUCN (2021), Protected Planet: The World Database on Protected Areas (WDPA) and World Database on Other Effective Area-based Conservation Measures (WD-OECM) [Online], August 2021, Cambridge, UK: UNEP-WCMC and IUCN. Available at: www.protectedplanet.net.
    - 2. Unzipp and rename the folder containing the .shp file to `protected_areas`
    - 3. Important! Don't delete the other files which come with the .shp file. They are required to build the shape.
    - 4. Move the file in such a way that the above path is given
    - 5. Activate the environment of environment-max.yaml
    - 6. Ready to run the script

Tip: The output file `natura.tiff` contains now the 100x100m rasters of protective areas. This operation can make the filesize of that TIFF quite large and leads to problems when trying to open. QGIS, an open source tool helps exploring the file.
"""
import logging
import os

import atlite
import geopandas as gpd
import rasterio as rio
from _helpers import configure_logging
from rasterio.features import geometry_mask
from rasterio.warp import transform_bounds

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)

CUTOUT_CRS = "EPSG:4326"


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


def determine_cutout_xXyY(cutout_name, out_logging):
    if out_logging:
        _logger.info("Stage 1/5: Determine cutout boundaries")
    cutout = atlite.Cutout(cutout_name)
    assert cutout.crs == CUTOUT_CRS
    x, X, y, Y = cutout.extent
    dx, dy = cutout.dx, cutout.dy
    return [x - dx / 2.0, X + dx / 2.0, y - dy / 2.0, Y + dy / 2.0]


def get_transform_and_shape(bounds, res, out_logging):
    if out_logging:
        _logger.info("Stage 2/5: Get transform and shape")
    left, bottom = [(b // res) * res for b in bounds[:2]]
    right, top = [(b // res + 1) * res for b in bounds[2:]]
    shape = int((top - bottom) // res), int((right - left) / res)
    transform = rio.Affine(res, 0, left, 0, -res, top)
    return transform, shape


def unify_protected_shape_areas(inputs, crs, out_logging):
    """
    Iterates thorugh all snakemake rule inputs and unifies shapefiles (.shp) only.

    The input is given in the Snakefile and shapefiles are given by .shp


    Returns
    -------
    unified_shape : GeoDataFrame with a unified "multishape"

    """
    import pandas as pd
    from shapely.ops import unary_union

    if out_logging:
        _logger.info("Stage 3/5: Unify protected shape area.")

    # Read only .shp snakemake inputs
    shp_files = [string for string in inputs if ".shp" in string]
    assert len(shp_files) != 0, "no input shapefiles given"
    # Create one geodataframe with all geometries, of all .shp files
    if out_logging:
        _logger.info(
            "Stage 3/5: Unify protected shape area. Step 1: Create one geodataframe with all shapes"
        )
    for i in shp_files:
        shape = gpd.GeoDataFrame(
            pd.concat([gpd.read_file(i) for i in shp_files])
        ).to_crs(crs)

    # Removes shapely geometry with null values. Returns geoseries.
    shape = shape["geometry"][shape["geometry"].is_valid]

    # Create Geodataframe with crs
    shape = gpd.GeoDataFrame(shape, crs=crs)
    shape = shape.rename(columns={0: "geometry"}).set_geometry("geometry")

    # Unary_union makes out of i.e. 1000 shapes -> 1 unified shape
    if out_logging:
        _logger.info("Stage 3/5: Unify protected shape area. Step 2: Unify all shapes")
    unified_shape_file = unary_union(shape["geometry"])
    if out_logging:
        _logger.info(
            "Stage 3/5: Unify protected shape area. Step 3: Set geometry of unified shape"
        )
    unified_shape = gpd.GeoDataFrame(geometry=[unified_shape_file], crs=crs)

    return unified_shape


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_natura_raster")
    configure_logging(snakemake)

    # get crs
    area_crs = snakemake.config["crs"]["area_crs"]

    out_logging = True
    inputs = snakemake.input
    cutouts = inputs.cutouts
    shapefiles = get_fileshapes(inputs)
    xs, Xs, ys, Ys = zip(
        *(determine_cutout_xXyY(cutout, out_logging=out_logging) for cutout in cutouts)
    )
    bounds = transform_bounds(CUTOUT_CRS, area_crs, min(xs), min(ys), max(Xs), max(Ys))
    transform, out_shape = get_transform_and_shape(
        bounds, res=100, out_logging=out_logging
    )
    # adjusted boundaries
    shapes = unify_protected_shape_areas(shapefiles, area_crs, out_logging=out_logging)

    if out_logging:
        _logger.info("Stage 4/5: Mask geometry")
    raster = ~geometry_mask(shapes.geometry, out_shape[::-1], transform)
    raster = raster.astype(rio.uint8)

    if out_logging:
        _logger.info("Stage 5/5: Export as .tiff")
    with rio.open(
        snakemake.output[0],
        "w",
        driver="GTiff",
        dtype=rio.uint8,
        count=1,
        transform=transform,
        crs=area_crs,
        compress="lzw",
        width=raster.shape[1],
        height=raster.shape[0],
    ) as dst:
        dst.write(raster, indexes=1)
