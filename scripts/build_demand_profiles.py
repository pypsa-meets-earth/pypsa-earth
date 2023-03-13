# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2023- PyPSA-Earth Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Creates electric demand profile csv

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

- ``networks/base.nc``: confer :ref:`base`, a base PyPSA Network 
- ``resources/bus_regions/regions_onshore.geojson``: confer :ref:`bus regions` 
- ``load_data_paths``: paths to load profiles, e.g. hourly country load profiles produced by GEGIS
- ``resources/shapes/gadm_shapes.geojson``: confer :ref:`shapes`, file containing the gadm shapes

Outputs
-------

- ``resources/demand_profiles.csv``: the content of the file is the electric demand profile associated to each bus. The file has the snapshots as rows and the buses of the network as columns. 

Description
-----------

The rule :mod:`build_demand` creates load demand profiles in correspondance of the buses of the network.
It creates the load paths for GEGIS outputs by combining the input parameters of the countries, weather year, prediction year, and SSP scenario.
Then with a function that takes in the PyPSA network "base.nc", region and gadm shape data, the countries of interest, a scale factor, and the snapshots,
it returns a csv file called "demand_profiles.csv", that allocates the load to the buses of the network according to GDP and population.

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


def get_load_paths_gegis(ssp_parentfolder, config):
    """
    Create load paths for GEGIS outputs

    The paths are created automatically according to included country,
    weather year, prediction year and ssp scenario

    Example
    -------
    ["/data/ssp2-2.6/2030/era5_2013/Africa.nc", "/data/ssp2-2.6/2030/era5_2013/Africa.nc"]
    """
    countries = config.get("countries")
    region_load = getContinent(countries)
    weather_year = config.get("load_options")["weather_year"]
    prediction_year = config.get("load_options")["prediction_year"]
    ssp = config.get("load_options")["ssp"]

    load_paths = []
    for continent in region_load:
        load_path = os.path.join(
            ssp_parentfolder,
            str(ssp),
            str(prediction_year),
            "era5_" + str(weather_year),
            str(continent) + ".nc",
        )
        load_paths.append(load_path)

    return load_paths


def build_demand_profiles(
    n,
    load_paths,
    regions,
    admin_shapes,
    countries,
    scale,
    start_date,
    end_date,
    out_path,
):
    """
    Create csv file of electric demand time series.

    Parameters
    ----------
    n : pypsa network
    load_paths: paths of the load files
    regions : .geojson
        Contains bus_id of low voltage substations and
        bus region shapes (voronoi cells)
    admin_shapes : .geojson
        contains subregional gdp, population and shape data
    countries : list
        List of countries that is config input
    scale : float
        The scale factor is multiplied with the load (1.3 = 30% more load)
    start_date: parameter
        The start_date is the first hour of the first day of the snapshots
    end_date: parameter
        The end_date is the last hour of the last day of the snapshots

    Returns
    -------
    demand_profiles.csv : csv file containing the electric demand time series
    """
    substation_lv_i = n.buses.index[n.buses["substation_lv"]]
    regions = gpd.read_file(regions).set_index("name").reindex(substation_lv_i)
    load_paths = load_paths
    # Merge load .nc files: https://stackoverflow.com/questions/47226429/join-merge-multiple-netcdf-files-using-xarray
    gegis_load = xr.open_mfdataset(load_paths, combine="nested")
    gegis_load = gegis_load.to_dataframe().reset_index().set_index("time")
    # filter load for analysed countries
    gegis_load = gegis_load.loc[gegis_load.region_code.isin(countries)]
    logger.info(f"Load data scaled with scaling factor {scale}.")
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

    demand_profiles = pd.concat(
        [
            upsample(cntry, group)
            for cntry, group in regions.geometry.groupby(regions.country)
        ],
        axis=1,
    )

    start_date = pd.to_datetime(start_date)
    end_date = pd.to_datetime(end_date) - pd.Timedelta(hours=1)
    demand_profiles = demand_profiles.loc[start_date:end_date]
    demand_profiles.to_csv(out_path, header=True)

    logger.info(f"Demand_profiles csv file created for the corresponding snapshots.")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_demand_profiles")
        sets_path_to_root("pypsa-earth")
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.base_network)

    # Snakemake imports:
    regions = snakemake.input.regions
    load_paths = snakemake.input["load"]
    countries = snakemake.config["countries"]
    admin_shapes = snakemake.input.gadm_shapes
    scale = snakemake.config["load_options"]["scale"]
    start_date = snakemake.config["snapshots"]["start"]
    end_date = snakemake.config["snapshots"]["end"]
    out_path = snakemake.output[0]

    build_demand_profiles(
        n,
        load_paths,
        regions,
        admin_shapes,
        countries,
        scale,
        start_date,
        end_date,
        out_path,
    )
