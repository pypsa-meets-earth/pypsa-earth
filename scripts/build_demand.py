# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors, 2021 PyPSA-Africa Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later
# coding: utf-8
"""
Adds electrical generators, load and existing hydro storage units to a base network.

Relevant Settings
-----------------

.. code:: yaml

    costs:
        year:
        version:
        rooftop_share:
        USD2013_to_EUR2013:
        dicountrate:
        emission_prices:

    electricity:
        max_hours:
        marginal_cost:
        capital_cost:
        conventional_carriers:
        co2limit:
        extendable_carriers:
        include_renewable_capacities_from_OPSD:
        estimate_renewable_capacities_from_capacity_stats:

    load:
        scale:
        ssp:
        weather_year:
        prediction_year:
        region_load:


    renewable:
        hydro:
            carriers:
            hydro_max_hours:
            hydro_capital_cost:

    lines:
        length_factor:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at :ref:`costs_cf`,
    :ref:`electricity_cf`, :ref:`load_cf`, :ref:`renewable_cf`, :ref:`lines_cf`

Inputs
------

- ``resources/costs.csv``: The database of cost assumptions for all included technologies for specific years from various sources; e.g. discount rate, lifetime, investment (CAPEX), fixed operation and maintenance (FOM), variable operation and maintenance (VOM), fuel costs, efficiency, carbon-dioxide intensity.
- ``data/bundle/hydro_capacities.csv``: Hydropower plant store/discharge power capacities, energy storage capacity, and average hourly inflow by country.  Not currently used!

    .. image:: ../img/hydrocapacities.png
        :scale: 34 %

- ``data/geth2015_hydro_capacities.csv``: alternative to capacities above; not currently used!
- ``resources/ssp2-2.6/2030/era5_2013/Africa.nc`` Hourly country load profiles produced by GEGIS
- ``resources/regions_onshore.geojson``: confer :ref:`busregions`
- ``resources/gadm_shapes.geojson``: confer :ref:`shapes`
- ``resources/powerplants.csv``: confer :ref:`powerplants`
- ``resources/profile_{}.nc``: all technologies in ``config["renewables"].keys()``, confer :ref:`renewableprofiles`.
- ``networks/base.nc``: confer :ref:`base`

Outputs
-------

- ``networks/elec.nc``:

    .. image:: ../img/elec.png
            :scale: 33 %

Description
-----------

The rule :mod:`add_electricity` ties all the different data inputs from the preceding rules together into a detailed PyPSA network that is stored in ``networks/elec.nc``. It includes:

- today's transmission topology and transfer capacities (in future, optionally including lines which are under construction according to the config settings ``lines: under_construction`` and ``links: under_construction``),
- today's thermal and hydro power generation capacities (for the technologies listed in the config setting ``electricity: conventional_carriers``), and
- today's load time-series (upsampled in a top-down approach according to population and gross domestic product)

It further adds extendable ``generators`` with **zero** capacity for

- photovoltaic, onshore and AC- as well as DC-connected offshore wind installations with today's locational, hourly wind and solar capacity factors (but **no** current capacities),
- additional open- and combined-cycle gas turbines (if ``OCGT`` and/or ``CCGT`` is listed in the config setting ``electricity: extendable_carriers``)
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

#idx = pd.IndexSlice

logger = logging.getLogger(__name__)


def normed(s):
    return s / s.sum()

# def get_load_paths_gegis(ssp_parentfolder, config):
#     """
#     Creates load paths for the GEGIS outputs

#     The paths are created automatically according to included country,
#     weather year, prediction year and ssp scenario

#     Example
#     -------
#     ["/data/ssp2-2.6/2030/era5_2013/Africa.nc", "/data/ssp2-2.6/2030/era5_2013/Africa.nc"]
#     """
#     countries = config.get("countries")
#     region_load = getContinent(countries)
#     weather_year = config.get("load_options")["weather_year"]
#     prediction_year = config.get("load_options")["prediction_year"]
#     ssp = config.get("load_options")["ssp"]

#     load_paths = []
#     for continent in region_load:
#         load_path = os.path.join(
#             ssp_parentfolder,
#             str(ssp),
#             str(prediction_year),
#             "era5_" + str(weather_year),
#             str(continent) + ".nc",
#         )
#         load_paths.append(load_path)

#     return load_paths


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
        snakemake = mock_snakemake("add_electricity")
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
