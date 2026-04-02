# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Creates Voronoi shapes for each bus representing both onshore and offshore
regions.

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

    .. image:: /img/regions_onshore.png
        :width: 33 %

- ``resources/regions_offshore.geojson``:

    .. image:: /img/regions_offshore.png
        :width: 33 %

Description
-----------
"""
import os
import warnings

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
from _helpers import REGION_COLS, configure_logging, create_logger, nearest_shape
from shapely.geometry import Polygon

logger = create_logger(__name__)


def voronoi(
    points: pd.DataFrame,
    outline: Polygon,
    geo_crs: str = "EPSG:4326",
) -> gpd.GeoSeries:
    """
    Create Voronoi polygons from a set of points within an outline.

    Parameters
    ----------
    points : pd.DataFrame
         DataFrame containing the coordinates of the points with columns ["x", "y"] and index
    outline : Polygon
        Shapely Polygon defining the outline within which to compute the Voronoi partition.
    geo_crs : str
        CRS used for geographic projection, passed to GeoPandas (e.g. "EPSG:4326")

    Returns
    -------
    gpd.GeoSeries
        GeoSeries of Voronoi polygons corresponding to each point in `points`, clipped to the `outline` polygon.
    """

    pts = gpd.GeoSeries(
        gpd.points_from_xy(points.x, points.y),
        index=points.index,
        crs=geo_crs,
    )
    voronoi = pts.voronoi_polygons(extend_to=outline).clip(outline)

    # can be removed with shapely 2.1 where order is preserved
    # https://github.com/shapely/shapely/issues/2020
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        pts = gpd.GeoDataFrame(geometry=pts)
        voronoi = gpd.GeoDataFrame(geometry=voronoi)
        joined = gpd.sjoin_nearest(pts, voronoi, how="right")

    gdf = joined.dissolve(by=points.index.name).reindex(points.index).squeeze()

    return gdf


def get_gadm_shape(
    onshore_buses: pd.DataFrame,
    gadm_shapes: gpd.GeoDataFrame,
    geo_crs: str = "EPSG:4326",
    metric_crs: str = "EPSG:3857",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Get the GADM shape for each bus by finding the nearest GADM shape to each bus.

    Parameters
    ----------
    onshore_buses: pd.DataFrame
        DataFrame containing the onshore buses with columns ["x", "y"]
    gadm_shapes: gpd.GeoDataFrame
        GeoDataFrame containing the GADM shapes with a geometry column
    geo_crs : str
        CRS used for geographic projection, passed to GeoPandas (e.g. "EPSG:4326")
    metric_crs : str
        CRS used for distance projection, passed to GeoPandas (e.g. "EPSG:3857")

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        A tuple containing two arrays:
        - An array of geometries corresponding to the GADM shapes for each bus.
        - An array of indices corresponding to the GADM shape IDs for each bus.
    """
    geo_regions = gpd.GeoDataFrame(
        onshore_buses[["x", "y"]],
        geometry=gpd.points_from_xy(onshore_buses["x"], onshore_buses["y"]),
        crs=geo_crs,
    ).to_crs(metric_crs)

    join_geos = gpd.sjoin_nearest(
        geo_regions, gadm_shapes.to_crs(metric_crs), how="left"
    )

    # when duplicates, keep only the first entry
    join_geos = join_geos[~join_geos.index.duplicated()]

    gadm_sel = gadm_shapes.loc[join_geos[gadm_shapes.index.name].values]

    return gadm_sel.geometry.values, gadm_sel.index.values


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_bus_regions")

    configure_logging(snakemake)

    inputs = snakemake.input
    country_shapes_fn = inputs.get("subregion_shapes") or inputs.country_shapes
    offshore_shapes_fn = inputs.get("subregion_offshore") or inputs.offshore_shapes
    countries = snakemake.params.countries
    geo_crs = snakemake.params.crs["geo_crs"]
    area_crs = snakemake.params.crs["area_crs"]
    metric_crs = snakemake.params.crs["distance_crs"]

    n = pypsa.Network(inputs.base_network)

    country_shapes = gpd.read_file(country_shapes_fn).set_index("name")["geometry"]

    offshore_shapes = gpd.read_file(offshore_shapes_fn)

    offshore_shapes = offshore_shapes.reindex(columns=REGION_COLS).set_index("name")[
        "geometry"
    ]

    # Option for subregion
    subregion_shapes = snakemake.input.get("subregion_shapes")
    if subregion_shapes:
        crs = {"geo_crs": geo_crs, "distance_crs": metric_crs}
        tolerance = snakemake.config.get("subregion", {}).get("tolerance", 100)
        n = nearest_shape(n, country_shapes_fn, crs, tolerance=tolerance)

        countries = list(country_shapes.index)

    gadm_shapes = gpd.read_file(inputs.gadm_shapes).set_index("GADM_ID")

    onshore_regions = []
    offshore_regions = []

    for country in countries:
        c_b = n.buses.country == country
        if n.buses.loc[c_b & n.buses.substation_lv, ["x", "y"]].empty:
            logger.warning(f"No low voltage buses found for {country}!")
            continue

        onshore_shape = country_shapes[country]
        onshore_locs = n.buses.loc[c_b & n.buses.substation_lv, ["x", "y"]]
        gadm_country = gadm_shapes[gadm_shapes.country == country]
        if snakemake.params.alternative_clustering:
            onshore_geometry, shape_id = get_gadm_shape(
                onshore_locs,
                gadm_country,
                geo_crs,
                metric_crs,
            )
        else:
            onshore_geometry = voronoi(onshore_locs, onshore_shape)
            shape_id = 0  # Not used

        temp_region = gpd.GeoDataFrame(
            {
                "name": onshore_locs.index,
                "x": onshore_locs["x"],
                "y": onshore_locs["y"],
                "geometry": onshore_geometry,
                "country": country,
                "shape_id": shape_id,
            },
            crs=geo_crs,
        )
        temp_region = temp_region[
            temp_region.geometry.is_valid & ~temp_region.geometry.is_empty
        ]
        onshore_regions.append(temp_region)

        # These two logging could be commented out
        if country not in offshore_shapes.index:
            logger.warning(f"No off-shore shapes for {country}")
            continue

        offshore_shape = offshore_shapes[country]

        if n.buses.loc[c_b & n.buses.substation_off, ["x", "y"]].empty:
            logger.warning(f"No off-shore substations found for {country}")
            continue
        else:
            offshore_locs = n.buses.loc[c_b & n.buses.substation_off, ["x", "y"]]
            shape_id = 0  # Not used
            offshore_geometry = voronoi(offshore_locs, offshore_shape)
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
            offshore_regions_c = offshore_regions_c[
                offshore_regions_c.geometry.is_valid
                & ~offshore_regions_c.geometry.is_empty
            ]
            offshore_regions.append(offshore_regions_c)

    # create geodataframe and remove nan shapes
    onshore_regions = gpd.GeoDataFrame(
        pd.concat(onshore_regions, ignore_index=True),
        crs=country_shapes.crs,
    ).dropna(axis="index", subset=["geometry"])

    if snakemake.params.alternative_clustering:
        # determine isolated buses
        n.determine_network_topology()
        non_isolated_buses = n.buses.duplicated(subset=["sub_network"], keep=False)
        isolated_buses = n.buses[~non_isolated_buses].index
        non_isolated_regions = onshore_regions[
            ~onshore_regions.name.isin(isolated_buses)
        ]
        isolated_regions = onshore_regions[onshore_regions.name.isin(isolated_buses)]

        # Combine regions while prioritizing non-isolated ones
        onshore_regions = pd.concat(
            [non_isolated_regions, isolated_regions]
        ).drop_duplicates("shape_id", keep="first")

        if len(onshore_regions) < len(gadm_country):
            logger.warning(
                f"The number of remaining of buses are less than the number of administrative clusters suggested!"
            )

    if subregion_shapes:
        logger.info("Deactivate subregion classificaition")
        original_shapes = snakemake.input.original_shapes
        n = nearest_shape(n, original_shapes, crs, tolerance=tolerance)

        onshore_regions["country"] = onshore_regions.name.map(n.buses.country)
        if offshore_regions:
            for offshore_region in offshore_regions:
                offshore_region["country"] = offshore_region.name.map(n.buses.country)

    onshore_regions = pd.concat([onshore_regions], ignore_index=True).to_file(
        snakemake.output.regions_onshore
    )

    if offshore_regions:
        # if a offshore_regions exists execute below
        pd.concat(offshore_regions, ignore_index=True).to_file(
            snakemake.output.regions_offshore
        )
    else:
        # if no offshore_regions exist save an empty offshore_shape
        offshore_shapes.to_frame().to_file(snakemake.output.regions_offshore)
