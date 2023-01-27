# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors, 2021 PyPSA-Africa Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later
# coding: utf-8
"""
Adds load to a network having generators and existing hydro storage units.

Relevant Settings
-----------------

.. code:: yaml

    load:
        scale:
        ssp:
        weather_year:
        prediction_year:
        region_load:


Inputs
------

- ``networks/elec.nc``: a network created in the script add_electricity with generators and hydro storage units,
- ``resources/bus_regions/regions_onshore.geojson``: confer :ref:`busregions`,
- ``load_data_paths``, 
- ``resources/shapes/gadm_shapes.geojson``: confer :ref:`shapes`,
- ``networks/elec.nc``: confer :ref:`elec`

Outputs
-------

- ``networks/elec_1.nc``:

    .. image:: ../img/elec_1.png
            :scale: 33 %

Description
-----------

The rule :mod:`build_demand` attaches the load to the network stored in networks/elec.nc and stores the results in ``networks/elec_1.nc``. It includes:

- today's load time-series (upsampled in a top-down approach according to population and gross domestic product)

"""


import logging
import os

import geopandas as gpd
import numpy as np
import pandas as pd
import powerplantmatching as pm
import pypsa
import xarray as xr
from _helpers import configure_logging, getContinent, update_p_nom_max
from shapely.validation import make_valid
from vresutils import transfer as vtransfer


logger = logging.getLogger(__name__)


def normed(s):
    return s / s.sum()

def attach_load(
    n,
    load_paths,
    regions,
    admin_shapes,
    countries,
    scale,
):
    """
    Add load to the network and distributes them according GDP and population.

    Parameters
    ----------
    n : pypsa network
    regions : .geojson
        Contains bus_id of low voltage substations and
        bus region shapes (voronoi cells)
    load_paths: paths of the load files
    admin_shapes : .geojson
        contains subregional gdp, population and shape data
    countries : list
        List of countries that is config input
    scale : float
        The scale factor is multiplied with the load (1.3 = 30% more load)

    Returns
    -------
    n : pypsa network
        Now attached with load time series
    """
    substation_lv_i = n.buses.index[n.buses["substation_lv"]]
    regions = gpd.read_file(regions).set_index("name").reindex(substation_lv_i)

    load_paths = load_paths
    # Merge load .nc files: https://stackoverflow.com/questions/47226429/join-merge-multiple-netcdf-files-using-xarray
    gegis_load = xr.open_mfdataset(load_paths, combine="nested")
    gegis_load = gegis_load.to_dataframe().reset_index().set_index("time")
    # filter load for analysed countries
    gegis_load = gegis_load.loc[gegis_load.region_code.isin(countries)]
    logger.info(f"Load data scaled with scalling factor {scale}.")
    gegis_load["Electricity demand"] *= scale
    shapes = gpd.read_file(admin_shapes).set_index("GADM_ID")
    shapes["geometry"] = shapes["geometry"].apply(lambda x: make_valid(x))

    def upsample(cntry, group):
        """
        Distributes load in country according to population and gdp
        """
        l = gegis_load.loc[gegis_load.region_code == cntry]["Electricity demand"]
        if len(group) == 1:
            return pd.DataFrame({group.index[0]: l})
        else:
            shapes_cntry = shapes.loc[shapes.country == cntry]
            transfer = vtransfer.Shapes2Shapes(
                group, shapes_cntry.geometry, normed=False
            ).T.tocsr()
            gdp_n = pd.Series(
                transfer.dot(shapes_cntry["gdp"].fillna(1.0).values), index=group.index
            )
            pop_n = pd.Series(
                transfer.dot(shapes_cntry["pop"].fillna(1.0).values), index=group.index
            )

            # relative factors 0.6 and 0.4 have been determined from a linear
            # regression on the country to EU continent load data
            # (refer to vresutils.load._upsampling_weights)
            # TODO: require adjustment for Africa
            factors = normed(0.6 * normed(gdp_n) + 0.4 * normed(pop_n))
            return pd.DataFrame(
                factors.values * l.values[:, np.newaxis],
                index=l.index,
                columns=factors.index,
            )

    load = pd.concat(
        [
            upsample(cntry, group)
            for cntry, group in regions.geometry.groupby(regions.country)
        ],
        axis=1,
    )

    n.madd("Load", substation_lv_i, bus=substation_lv_i, p_set=load)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_demand")
        sets_path_to_root("pypsa-earth")
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.elec_network)

    # Snakemake imports:
    regions = snakemake.input.regions

    load_paths = snakemake.input["load"]
    countries = snakemake.config["countries"]
    admin_shapes = snakemake.input.gadm_shapes
    scale = snakemake.config["load_options"]["scale"]

    attach_load(n, load_paths, regions, admin_shapes, countries, scale)

    n.meta = snakemake.config
    n.export_to_netcdf(snakemake.output[0])
