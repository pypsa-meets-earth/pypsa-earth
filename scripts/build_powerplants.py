# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Retrieves conventional powerplant capacities and locations from `powerplantmatching <https://github.com/FRESNA/powerplantmatching>`_, assigns these to buses and creates a ``.csv`` file. It is possible to amend the powerplant database with custom entries provided in ``data/custom_powerplants.csv``.

Relevant Settings
-----------------

.. code:: yaml

    electricity:
      powerplants_filter:
      custom_powerplants:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`electricity`

Inputs
------

- ``networks/base.nc``: confer :ref:`base`.
- ``data/custom_powerplants.csv``: custom powerplants in the same format as `powerplantmatching <https://github.com/FRESNA/powerplantmatching>`_ provides or as OSM extractor generates

Outputs
-------

- ``resource/powerplants.csv``: A list of conventional power plants (i.e. neither wind nor solar) with fields for name, fuel type, technology, country, capacity in MW, duration, commissioning year, retrofit year, latitude, longitude, and dam information as documented in the `powerplantmatching README <https://github.com/FRESNA/powerplantmatching/blob/master/README.md>`_; additionally it includes information on the closest substation/bus in ``networks/base.nc``.

    .. image:: /img/powerplantmatching.png
        :width: 30 %

    **Source:** `powerplantmatching on GitHub <https://github.com/FRESNA/powerplantmatching>`_

Description
-----------

The configuration options ``electricity: powerplants_filter`` and ``electricity: custom_powerplants`` can be used to control whether data should be retrieved from the original powerplants database or from custom amendmends. These specify `pandas.query <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.query.html>`_ commands.

1. Adding all powerplants from custom:

    .. code:: yaml

        powerplants_filter: false
        custom_powerplants: true

2. Replacing powerplants in e.g. Germany by custom data:

    .. code:: yaml

        powerplants_filter: Country not in ['Germany']
        custom_powerplants: true

    or

    .. code:: yaml

        powerplants_filter: Country not in ['Germany']
        custom_powerplants: Country in ['Germany']

3. Adding additional built year constraints:

    .. code:: yaml

        powerplants_filter: Country not in ['Germany'] and YearCommissioned <= 2015
        custom_powerplants: YearCommissioned <= 2015

Format required for the custom_powerplants.csv should be similar to the powerplantmatching format with some additional considerations:

Columns required: [id, Name, Fueltype, Technology, Set, Country, Capacity, Efficiency, DateIn, DateRetrofit, DateOut, lat, lon, Duration, Volume_Mm3, DamHeight_m, StorageCapacity_MWh, EIC, projectID]

Tagging considerations for columns in the file:

- FuelType: 'Natural Gas' has to be tagged either as 'OCGT', 'CCGT'
- Technology: 'Reservoir' has to be set as 'ror' if hydro powerplants are to be considered as 'Generators' and not 'StorageUnits'
- Country:  Country name has to be defined with its alpha2 code ('NG' for Nigeria,'BO' for Bolivia, 'FR' for France, etc.

The following assumptions were done to map custom OSM-extracted power plants with powerplantmatching format.

1. The benchmark PPM keys values were taken as follows:
        'Fueltype': ['Hydro', 'Hard Coal', 'Natural Gas', 'Lignite', 'Nuclear', 'Oil', 'Bioenergy'
            'Wind', 'Geothermal', 'Solar', 'Waste', 'Other']

        'Technology': ['Reservoir', 'Pumped Storage', 'Run-Of-River', 'Steam Turbine', 'CCGT', 'OCGT'
            'Pv', 'CCGT, Thermal', 'Offshore', 'Storage Technologies']

        'Set': ['Store', 'PP', 'CHP']

2. OSM-extracted features were mapped into PPM ones using a (quite arbitrary) set of rules:
        'coal': 'Hard Coal'
        'wind_turbine': 'Onshore',
        'horizontal_axis' : 'Onshore',
        'vertical_axis' : 'Offhore',
        'nuclear': 'Steam Turbine'
3. All hydro OSM-extracted objects were interpreted as generation technologies, although ["Run-Of-River", "Pumped Storage", "Reservoir"] in PPM can belong to 'Storage Technologies', too.
4. OSM extraction was supposed to be ignoring non-generation features like CHP and Natural Gas storage (in contrast to PPM).
"""

import geopandas as gpd
import numpy as np
import pandas as pd
import powerplantmatching as pm
import pypsa
import yaml
from _helpers import (
    change_to_script_dir,
    configure_logging,
    create_logger,
    get_current_directory_path,
    get_path,
    mock_snakemake,
    read_csv_nafix,
    two_digits_2_name_country,
)
from scipy.spatial import cKDTree as KDTree
from shapely.geometry import Point

logger = create_logger(__name__)


def add_power_plants(
    custom_power_plants_file_path,
    power_plants_config_dict,
    ppl_assignment_strategy,
    countries_names_list,
):

    if ppl_assignment_strategy == "replace":
        # use only the data from custom_powerplants.csv
        custom_power_plants = read_csv_nafix(
            custom_power_plants_file_path, index_col=0, dtype={"bus": "str"}
        )
        return custom_power_plants
    elif ppl_assignment_strategy == "merge":
        # merge the data from powerplantmatching and custom_powerplants.csv
        ppl_ppm = (
            pm.powerplants(
                from_url=False, update=True, config_update=power_plants_config_dict
            )
            .powerplant.fill_missing_decommissioning_years()
            .query(
                'Fueltype not in ["Solar", "Wind"] and Country in @countries_names_list'
            )
            .powerplant.convert_country_to_alpha2()
            .pipe(replace_natural_gas_technology)
        )
        ppl_cpp = read_csv_nafix(
            custom_power_plants_file_path, index_col=0, dtype={"bus": "str"}
        )
        power_plants = pd.concat(
            [ppl_ppm, ppl_cpp], sort=False, ignore_index=True, verify_integrity=True
        )
        return power_plants
    elif (
        ppl_assignment_strategy not in ["merge", "replace"]
        or ppl_assignment_strategy is None
    ):
        # use only the data from powerplantsmatching
        power_plants = (
            pm.powerplants(
                from_url=False, update=True, config_update=power_plants_config_dict
            )
            .powerplant.fill_missing_decommissioning_years()
            .query(
                'Fueltype not in ["Solar", "Wind"] and Country in @countries_names_list'
            )
            .powerplant.convert_country_to_alpha2()
            .pipe(replace_natural_gas_technology)
        )
        return power_plants
    else:
        raise Exception(
            "No power plants were built for custom_powerplants {}".format(
                ppl_assignment_strategy
            )
        )


def replace_natural_gas_technology(df: pd.DataFrame):
    """
    Maps and replaces gas technologies in the powerplants.csv onto model
    compliant carriers.
    """
    mapping = {
        "Steam Turbine": "CCGT",
        "Combustion Engine": "OCGT",
        "NG": "CCGT",
        "Ng": "CCGT",
        "NG/FO": "OCGT",
        "Ng/Fo": "OCGT",
        "NG/D": "OCGT",
        "LNG": "OCGT",
        "CCGT/D": "CCGT",
        "CCGT/FO": "CCGT",
        "LCCGT": "CCGT",
        "CCGT/Fo": "CCGT",
    }
    fueltype = df["Fueltype"] == "Natural Gas"
    df.loc[fueltype, "Technology"] = (
        df.loc[fueltype, "Technology"].replace(mapping).fillna("CCGT")
    )
    unique_tech_with_ng = df.loc[fueltype, "Technology"].unique()
    unknown_techs = np.setdiff1d(unique_tech_with_ng, ["CCGT", "OCGT"])
    if len(unknown_techs) > 0:
        df.Technology.where(
            fueltype,
            df["Technology"].map({t: "CCGT" for t in unknown_techs}),
            inplace=True,
        )
    df["Fueltype"] = np.where(fueltype, df["Technology"], df["Fueltype"])
    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        change_to_script_dir(__file__)
        snakemake = mock_snakemake("build_powerplants")

    configure_logging(snakemake)

    with open(snakemake.input.pm_config, "r") as f:
        powerplants_config = yaml.safe_load(f)

    filepath_osm_ppl = snakemake.input.osm_powerplants
    filepath_osm2pm_ppl = snakemake.output.powerplants_osm2pm
    powerplants_assignment_strategy = snakemake.params.custom_powerplants

    n = pypsa.Network(snakemake.input.base_network)
    countries_codes = n.buses.country.unique()
    countries_names = list(map(two_digits_2_name_country, countries_codes))

    powerplants_config["target_countries"] = countries_names

    if (
        "EXTERNAL_DATABASE"
        in powerplants_config["matching_sources"]
        + powerplants_config["fully_included_sources"]
    ):
        if "EXTERNAL_DATABASE" not in powerplants_config:
            logger.error(
                "Missing configuration EXTERNAL_DATABASE in powerplantmatching config yaml\n\t"
                "Please check file configs/powerplantmatching_config.yaml"
            )
        logger.info("Parsing OSM generator data to powerplantmatching format")
        powerplants_config["EXTERNAL_DATABASE"]["fn"] = get_path(
            get_current_directory_path(), filepath_osm2pm_ppl
        )
    else:
        # create an empty file
        with open(filepath_osm2pm_ppl, "w"):
            pass

    # specify the main query for filtering powerplants
    ppl_query = snakemake.params.powerplants_filter
    if isinstance(ppl_query, str):
        powerplants_config["main_query"] = ppl_query
    else:
        powerplants_config["main_query"] = ""

    ppl = add_power_plants(
        snakemake.input.custom_powerplants_file,
        powerplants_config,
        powerplants_assignment_strategy,
        countries_names,
    )

    countries_without_ppl = [
        c for c in countries_codes if c not in ppl.Country.unique()
    ]

    for c in countries_codes:
        substation_i = n.buses.query("substation_lv and country == @c").index
        kdtree = KDTree(n.buses.loc[substation_i, ["x", "y"]].values)
        ppl_i = ppl.query("Country == @c").index

        tree_i = kdtree.query(ppl.loc[ppl_i, ["lon", "lat"]].values)[1]
        ppl.loc[ppl_i, "bus"] = substation_i.append(pd.Index([np.nan]))[tree_i]

    if countries_without_ppl:
        logger.warning(f"No powerplants known in: {', '.join(countries_without_ppl)}")

    bus_null_b = ppl["bus"].isnull()
    if bus_null_b.any():
        logger.warning(f"Couldn't find close bus for {bus_null_b.sum()} powerplants")

    if snakemake.params.alternative_clustering:
        gadm_layer_id = snakemake.params.gadm_layer_id
        country_list = snakemake.params.countries
        geo_crs = snakemake.params.geo_crs

        gdf = gpd.read_file(snakemake.input.gadm_shapes)

        def locate_bus(coords, co):
            gdf_co = gdf[gdf["GADM_ID"].str.contains(co)]

            point = Point(coords["lon"], coords["lat"])

            try:
                return gdf_co[gdf_co.contains(point)][
                    "GADM_ID"
                ].item()  # filter gdf_co which contains point and returns the bus

            except ValueError:
                return gdf_co[
                    gdf_co.geometry == min(gdf_co.geometry, key=(point.distance))
                ][
                    "GADM_ID"
                ].item()  # looks for closest one shape=node
                # fixing https://github.com/pypsa-meets-earth/pypsa-earth/pull/670

        ppl["region_id"] = ppl[["lon", "lat", "Country"]].apply(
            lambda pp: locate_bus(pp[["lon", "lat"]], pp["Country"]), axis=1
        )

    ppl.to_csv(snakemake.output.powerplants)
