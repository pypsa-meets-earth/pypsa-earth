# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors, 2021 PyPSA-Africa Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Creates Voronoi shapes for each bus representing both onshore and offshore regions.

Relevant Settings
-----------------

.. code:: yaml

    countries:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`toplevel_cf`

Inputs
------

- ``resources/country_shapes.geojson``: confer :ref:`shapes`
- ``resources/offshore_shapes.geojson``: confer :ref:`shapes`
- ``networks/base.nc``: confer :ref:`base`

Outputs
-------

- ``resources/regions_onshore.geojson``:

    .. image:: ../img/regions_onshore.png
        :scale: 33 %

- ``resources/regions_offshore.geojson``:

    .. image:: ../img/regions_offshore.png
        :scale: 33 %

Description
-----------

"""
import logging
import os

import geopandas as gpd
import numpy
import pandas as pd
import pypsa
from _helpers import configure_logging, two_2_three_digits_country
from shapely.geometry import Point, Polygon
from vresutils.graph import voronoi_partition_pts

# from scripts.build_shapes import gadm

_logger = logging.getLogger(__name__)


def save_to_geojson(df, fn):
    # remove file if it exists
    if os.path.exists(fn):
        os.unlink(fn)
    if not isinstance(df, gpd.GeoDataFrame):
        df = gpd.GeoDataFrame(dict(geometry=df))

    # save file if the GeoDataFrame is non-empty
    if df.shape[0] > 0:
        df = df.reset_index()
        schema = {**gpd.io.file.infer_schema(df), "geometry": "Unknown"}
        df.to_file(fn, driver="GeoJSON", schema=schema)
    else:
        # create empty file to avoid issues with snakemake
        with open(fn, "w") as fp:
            pass


def custom_voronoi_partition_pts(points, outline, add_bounds_shape=True, multiplier=5):
    """
    Compute the polygons of a voronoi partition of `points` within the
    polygon `outline`

    Attributes
    ----------
    points : Nx2 - ndarray[dtype=float]
    outline : Polygon
    no_multipolygons : bool (default: False)
        If true, replace each MultiPolygon by its largest component

    Returns
    -------
    polygons : N - ndarray[dtype=Polygon|MultiPolygon]
    """

    import numpy as np
    from scipy.spatial import Voronoi
    from shapely.geometry import Point, Polygon

    points = np.asarray(points)

    polygons_arr = []

    if len(points) == 1:
        polygons_arr = [outline]
    else:

        xmin, ymin = np.amin(points, axis=0)
        xmax, ymax = np.amax(points, axis=0)

        if add_bounds_shape:
            # check bounds of the shape
            minx_o, miny_o, maxx_o, maxy_o = outline.boundary.bounds
            xmin = min(xmin, minx_o)
            ymin = min(ymin, miny_o)
            xmax = min(xmax, maxx_o)
            ymax = min(ymax, maxy_o)

        xspan = xmax - xmin
        yspan = ymax - ymin

        # to avoid any network positions outside all Voronoi cells, append
        # the corners of a rectangle framing these points
        vor = Voronoi(
            np.vstack(
                (
                    points,
                    [
                        [xmin - multiplier * xspan, ymin - multiplier * yspan],
                        [xmin - multiplier * xspan, ymax + multiplier * yspan],
                        [xmax + multiplier * xspan, ymin - multiplier * yspan],
                        [xmax + multiplier * xspan, ymax + multiplier * yspan],
                    ],
                )
            )
        )

        polygons_arr = np.empty((len(points),), "object")
        for i in range(len(points)):
            poly = Polygon(vor.vertices[vor.regions[vor.point_region[i]]])

            if not poly.is_valid:
                poly = poly.buffer(0)

            poly = poly.intersection(outline)

            polygons_arr[i] = poly

    return polygons_arr


def get_gadm_shape(onshore_locs, gadm_shapes):
    def locate_bus(coords):
        point = Point(Point(coords["x"], coords["y"]))
        gadm_shapes_country = gadm_shapes.filter(
            like=two_2_three_digits_country(country), axis=0
        )

        try:
            return gadm_shapes[gadm_shapes.contains(point)].item()
        except ValueError:
            return min(
                gadm_shapes_country, key=(point.distance)
            )  # TODO returns closest shape if the point was not inside one. Works well but will not catch an outlier bus.

    def get_id(coords):
        point = Point(Point(coords["x"], coords["y"]))
        gadm_shapes_country = gadm_shapes.filter(
            like=two_2_three_digits_country(country), axis=0
        )

        try:
            return gadm_shapes[gadm_shapes.contains(point)].index.item()
        except ValueError:
            # return 'not_found'
            return gadm_shapes_country[
                gadm_shapes_country.geometry
                == min(gadm_shapes_country, key=(point.distance))
            ].index.item()  # TODO returns closest shape if the point was not inside one.

    regions = onshore_locs[["x", "y"]].apply(locate_bus, axis=1)
    ids = onshore_locs[["x", "y"]].apply(get_id, axis=1)
    return regions.values, ids.values


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_bus_regions")
    configure_logging(snakemake)

    countries = snakemake.config["countries"]
    area_crs = snakemake.config["crs"]["area_crs"]

    n = pypsa.Network(snakemake.input.base_network)

    country_shapes = gpd.read_file(snakemake.input.country_shapes).set_index("name")[
        "geometry"
    ]
    offshore_shapes = gpd.read_file(snakemake.input.offshore_shapes).set_index("name")[
        "geometry"
    ]
    gadm_shapes = gpd.read_file(snakemake.input.gadm_shapes).set_index("GADM_ID")[
        "geometry"
    ]

    onshore_regions = []
    offshore_regions = []

    for country in countries:

        c_b = n.buses.country == country
        if n.buses.loc[c_b & n.buses.substation_lv, ["x", "y"]].empty:
            _logger.warning(f"No low voltage buses found for {country}!")
            continue

        onshore_shape = country_shapes[country]
        onshore_locs = n.buses.loc[c_b & n.buses.substation_lv, ["x", "y"]]
        if snakemake.config["cluster_options"]["alternative_clustering"]:
            onshore_geometry = get_gadm_shape(onshore_locs, gadm_shapes)[0]
            shape_id = get_gadm_shape(onshore_locs, gadm_shapes)[1]
        else:
            onshore_geometry = custom_voronoi_partition_pts(
                onshore_locs.values, onshore_shape
            )
            shape_id = 0  # Not used
        onshore_regions.append(
            gpd.GeoDataFrame(
                {
                    "name": onshore_locs.index,
                    "x": onshore_locs["x"],
                    "y": onshore_locs["y"],
                    "geometry": onshore_geometry,
                    "country": country,
                    "shape_id": shape_id,
                },
                crs=country_shapes.crs,
            )
        )

        # These two logging could be commented out
        if country not in offshore_shapes.index:
            _logger.warning(f"No off-shore shapes for {country}")
            continue

        offshore_shape = offshore_shapes[country]

        if n.buses.loc[c_b & n.buses.substation_off, ["x", "y"]].empty:
            _logger.warning(f"No off-shore substations found for {country}")
            continue
        else:
            offshore_locs = n.buses.loc[c_b & n.buses.substation_off, ["x", "y"]]
            shape_id = 0  # Not used
            offshore_geometry = custom_voronoi_partition_pts(
                offshore_locs.values, offshore_shape
            )
            offshore_regions_c = gpd.GeoDataFrame(
                {
                    "name": offshore_locs.index,
                    "x": offshore_locs["x"],
                    "y": offshore_locs["y"],
                    "geometry": offshore_geometry,
                    "country": country,
                    "shape_id": shape_id,
                },
                crs=country_shapes.crs,
            )
            offshore_regions_c = offshore_regions_c.loc[
                offshore_regions_c.to_crs(area_crs).area > 1e-2
            ]
            offshore_regions.append(offshore_regions_c)

    # create geodataframe and remove nan shapes
    onshore_regions = gpd.GeoDataFrame(
        pd.concat(onshore_regions, ignore_index=True),
        crs=country_shapes.crs,
    ).dropna(axis="index", subset=["geometry"])

    # remove empty polygons
    onshore_regions = onshore_regions[~onshore_regions.geometry.is_empty]

    save_to_geojson(
        onshore_regions,
        snakemake.output.regions_onshore,
    )
    if len(offshore_regions) != 0:
        # create geodataframe and remove nan shapes
        offshore_regions = gpd.GeoDataFrame(
            pd.concat(offshore_regions, ignore_index=True),
            crs=country_shapes.crs,
        ).dropna(axis="index", subset=["geometry"])

        # remove empty polygons
        offshore_regions = offshore_regions[~offshore_regions.geometry.is_empty]

    save_to_geojson(offshore_regions, snakemake.output.regions_offshore)
