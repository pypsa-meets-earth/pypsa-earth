# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
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
Steps to retrieve the protected area data (as apparently no API is given for the WDPA data):
    - 1. Download the WPDA Dataset: World Database on Protected Areas. UNEP-WCMC and IUCN (2021), Protected Planet: The World Database on Protected Areas (WDPA) and World Database on Other Effective Area-based Conservation Measures (WD-OECM) [Online], August 2021, Cambridge, UK: UNEP-WCMC and IUCN. Available at: www.protectedplanet.net.
    - 2. Unzipp and rename the folder containing the .shp file to `protected_areas`
    - 3. Important! Don't delete the other files which come with the .shp file. They are required to build the shape.
    - 4. Move the file in such a way that the above path is given
    - 5. Activate the environment of environment-max.yaml
    - 6. Ready to run the script

Tip: The output file `natura.tiff` contains now the 100x100m rasters of protective areas. This operation can make the filesize of that TIFF quite large and leads to problems when trying to open. QGIS, an open source tool helps exploring the file.
"""
import os
import logging
from _helpers import configure_logging

import atlite
import geopandas as gpd
import rasterio as rio
from rasterio.features import geometry_mask
from rasterio.warp import transform_bounds

logger = logging.getLogger(__name__)

# Requirement to set path to filepath for execution
os.chdir(os.path.dirname(os.path.abspath(__file__)))


def determine_cutout_xXyY(cutout_name):
    cutout = atlite.Cutout(cutout_name)
    assert cutout.crs.to_epsg() == 4326
    x, X, y, Y = cutout.extent
    dx, dy = cutout.dx, cutout.dy
    return [x - dx/2., X + dx/2., y - dy/2., Y + dy/2.]


def get_transform_and_shape(bounds, res):
    left, bottom = [(b // res)* res for b in bounds[:2]]
    right, top = [(b // res + 1) * res for b in bounds[2:]]
    shape = int((top - bottom) // res), int((right - left) / res)
    transform = rio.Affine(res, 0, left, 0, -res, top)
    return transform, shape


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_natura_raster')
    configure_logging(snakemake)


    cutouts = snakemake.input.cutouts
    xs, Xs, ys, Ys = zip(*(determine_cutout_xXyY(cutout) for cutout in cutouts))
    bounds = transform_bounds(4326, 3035, min(xs), min(ys), max(Xs), max(Ys))
    transform, out_shape = get_transform_and_shape(bounds, res=100)

    # adjusted boundaries
    shapes = gpd.read_file(snakemake.input.natura).to_crs(3035)
    raster = ~geometry_mask(shapes.geometry, out_shape[::-1], transform)
    raster = raster.astype(rio.uint8)

    with rio.open(snakemake.output[0], 'w', driver='GTiff', dtype=rio.uint8,
                  count=1, transform=transform, crs=3035, compress='lzw',
                  width=raster.shape[1], height=raster.shape[0]) as dst:
        dst.write(raster, indexes=1)

