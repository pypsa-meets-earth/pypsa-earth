# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
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
from _helpers import configure_logging

import pypsa
import os
import pandas as pd
import geopandas as gpd

from vresutils.graph import voronoi_partition_pts

logger = logging.getLogger(__name__)

# Requirement to set path to filepath for execution
# os.chdir(os.path.dirname(os.path.abspath(__file__)))


def save_to_geojson(s, fn):
    if os.path.exists(fn):
        os.unlink(fn)
    schema = {**gpd.io.file.infer_schema(s), 'geometry': 'Unknown'}
    s.to_file(fn, driver='GeoJSON', schema=schema)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_bus_regions')
    configure_logging(snakemake)

    countries = snakemake.config['countries']

    n = pypsa.Network(snakemake.input.base_network)

    country_shapes = gpd.read_file(snakemake.input.country_shapes).set_index('name')['geometry']
    offshore_shapes = gpd.read_file(snakemake.input.offshore_shapes).set_index('name')['geometry']

    # Issues in voronoi_creation due to overlapping/duplications will be removed with that function
    # First issues where observed for offshore shapes means that first country-onshore voronoi worked fine
    # TODO: Find out the root cause for this issues
    import shapely
    from shapely.validation import make_valid
    from shapely.geometry import shape, JOIN_STYLE
    from shapely.ops import unary_union
    # country_shapes = make_valid(country_shapes)
    # offshore_shapes = make_valid(offshore_shapes)

    onshore_regions = []
    offshore_regions = []

    for country in countries:

        c_b = n.buses.country == country
        # Check if lv_buses exist in onshape. TD has no buses!
        if n.buses.loc[c_b & n.buses.substation_lv, ["x", "y"]].empty: continue # TODO: Log warning! Something is going wrong.
        print(country)
        onshore_shape = make_valid(country_shapes[country])
        print(shapely.validation.explain_validity(onshore_shape), onshore_shape.area)
        onshore_locs = n.buses.loc[c_b & n.buses.substation_lv, ["x", "y"]]
        print(onshore_locs.values)
        onshore_regions.append(gpd.GeoDataFrame({
                'name': onshore_locs.index,
                'x': onshore_locs['x'],
                'y': onshore_locs['y'],
                'geometry': voronoi_partition_pts(onshore_locs.values, onshore_shape),
                'country': country
            }))

        if country not in offshore_shapes.index: continue
        if n.buses.loc[c_b & n.buses.substation_off, ["x", "y"]].empty: continue
        # TODO: fix issues 
        if country == "CD": continue # TODO: Remove error for Voronoi aggregation..
        if country == "KE": continue # TODO: Remove validity error..
        if country == "SD": continue # TODO: Remove validity error..

        offshore_shape = offshore_shapes[country]
        if offshore_shape.is_empty: continue

        offshore_shape = make_valid(offshore_shape) # Issue with CM reqired buffer
        print(offshore_shape.is_valid)
        print(offshore_shape.is_simple)
        print(shapely.validation.explain_validity(offshore_shape), offshore_shape.area)
        offshore_locs = n.buses.loc[c_b & n.buses.substation_off, ["x", "y"]]
        offshore_regions_c = gpd.GeoDataFrame({
                'name': offshore_locs.index,
                'x': offshore_locs['x'],
                'y': offshore_locs['y'],
                'geometry': voronoi_partition_pts(offshore_locs.values, offshore_shape),
                'country': country
            })
        offshore_regions_c = offshore_regions_c.loc[offshore_regions_c.area > 1e-2]
        offshore_regions.append(offshore_regions_c)

    save_to_geojson(pd.concat(onshore_regions, ignore_index=True), snakemake.output.regions_onshore)
    if len(offshore_regions) != 0: 
        offshore_regions= pd.concat(offshore_regions, ignore_index=True)
    save_to_geojson(offshore_regions, snakemake.output.regions_offshore)
