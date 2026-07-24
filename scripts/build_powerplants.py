# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Retrieves conventional powerplant capacities and locations from `powerplantmatching <https://github.com/FRESNA/powerplantmatching>`_, assigns these to buses and creates a ``.csv`` file. It is possible to amend or replace country-specific entries in the powerplant database with custom entries provided through one or more custom powerplant files.

Relevant Settings
-----------------

.. code:: yaml

    electricity:
      powerplants_filter:
      custom_powerplants:
        filepaths:
        method:

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

The configuration option ``electricity: powerplants_filter`` specifies a `pandas.query <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.query.html>`_ command applied to the original powerplantmatching database.

The ``electricity: custom_powerplants`` section specifies the custom files and how they are applied. ``method`` accepts ``false``, ``merge``, or ``replace``. It may be a single value applied to every filepath or a list matching ``filepaths``. The countries affected by ``replace`` are determined from the ``Country`` column of the corresponding file.

1. Using only powerplantmatching data:

    .. code:: yaml

        custom_powerplants:
          filepaths:
          - data/custom_powerplants.csv
          method: false

2. Adding all powerplants from a custom file:

    .. code:: yaml

        custom_powerplants:
          filepaths:
          - data/custom_powerplants.csv
          method: merge

3. Replacing powerplants for the countries contained in a custom file:

    .. code:: yaml

        custom_powerplants:
          filepaths:
          - data/custom_powerplants_US.csv
          method: replace

4. Applying different methods to different custom files:

    .. code:: yaml

        custom_powerplants:
          filepaths:
          - data/custom_powerplants_US.csv
          - data/custom_powerplants_BR.csv
          method:
          - replace
          - merge

Format required for the custom_powerplants.csv should be similar to the powerplantmatching format with some additional considerations:

Columns required: [id, Name, Fueltype, Technology, Set, Country, Capacity, Efficiency, DateIn, DateRetrofit, DateOut, lat, lon, Duration, Volume_Mm3, DamHeight_m, StorageCapacity_MWh, EIC, projectID]

Tagging considerations for columns in the file:

- FuelType: 'Natural Gas' has to be tagged either as 'OCGT', 'CCGT'
- Technology: 'Reservoir' has to be set as 'ror' if hydro powerplants are to be considered as 'Generators' and not 'StorageUnits'
- Country: Country name has to be defined with its alpha2 code ('NG' for Nigeria, 'BO' for Bolivia, 'FR' for France, etc.

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

import os

import geopandas as gpd
import numpy as np
import pandas as pd
import powerplantmatching as pm
import pypsa
import yaml
from _helpers import (
    configure_logging,
    create_logger,
    locate_bus,
    read_csv_nafix,
    to_csv_nafix,
    two_digits_2_name_country,
)
from scipy.spatial import cKDTree as KDTree
from shapely.geometry import Point

logger = create_logger(__name__)


def convert_osm_to_pm(filepath_ppl_osm, filepath_ppl_pm):
    if os.stat(filepath_ppl_osm).st_size == 0:
        return to_csv_nafix(pd.DataFrame(), filepath_ppl_pm, index=False)

    add_ppls = read_csv_nafix(filepath_ppl_osm, index_col=0, dtype={"bus": "str"})

    custom_ppls_coords = gpd.GeoSeries.from_wkt(add_ppls["geometry"])
    add_ppls = (
        add_ppls.rename(
            columns={
                "name": "Name",
                "tags.generator:source": "Fueltype",
                "tags.generator:type": "Technology",
                "tags.power": "Set",
                "power_output_MW": "Capacity",
            }
        )
        .replace(
            dict(
                Fueltype={
                    "nuclear": "Nuclear",
                    "wind": "Wind",
                    "hydro": "Hydro",
                    "tidal": "Other",
                    "wave": "Other",
                    "geothermal": "Geothermal",
                    "solar": "Solar",
                    # "Hard Coal" follows defaults of PPM
                    "coal": "Hard Coal",
                    "gas": "Natural Gas",
                    "biomass": "Bioenergy",
                    "biofuel": "Bioenergy",
                    "biogas": "Bioenergy",
                    "oil": "Oil",
                    "diesel": "Oil",
                    "gasoline": "Oil",
                    "waste": "Waste",
                    "osmotic": "Other",
                    "wave": "Other",
                    # approximation
                    # TODO: this shall be improved, one entry shall be Oil and the otherone gas
                    "gas;oil": "Oil",
                    "steam": "Natural Gas",
                    "waste_heat": "Other",
                },
                Technology={
                    "combined_cycle": "CCGT",
                    "gas_turbine": "OCGT",
                    "steam_turbine": "Steam Turbine",
                    "reciprocating_engine": "Combustion Engine",
                    # a very strong assumption
                    "wind_turbine": "Onshore",
                    "horizontal_axis": "Onshore",
                    "vertical_axis": "Offhore",
                    "solar_photovoltaic_panel": "Pv",
                },
                Set={"generator": "PP", "plant": "PP"},
            )
        )
        .assign(
            Country=lambda df: df.Country.map(two_digits_2_name_country),
            # Name=lambda df: "OSM_"
            # + df.Country.astype(str)
            # + "_"
            # + df.id.astype(str)
            # + "-"
            # + df.Name.astype(str),
            Efficiency="",
            Duration="",
            Volume_Mm3="",
            DamHeight_m="",
            StorageCapacity_MWh="",
            DateIn="",
            DateRetrofit="",
            DateMothball="",
            DateOut="",
            lat=custom_ppls_coords.y,
            lon=custom_ppls_coords.x,
            EIC=lambda df: df.id,
            projectID=lambda df: "OSM" + df.id.astype(str),
        )
        .dropna(subset=["Fueltype"])
    )

    # All Hydro objects can be interpreted by PPM as Storages, too
    # However, everything extracted from OSM seems to belong
    # to power plants with "tags.power" == "generator" only
    osm_ppm_df = pd.DataFrame(
        data={
            "osm_method": ["run-of-the-river", "water-pumped-storage", "water-storage"],
            "ppm_technology": ["Run-Of-River", "Pumped Storage", "Reservoir"],
        }
    )
    for i in osm_ppm_df.index:
        add_ppls.loc[
            add_ppls["tags.generator:method"] == osm_ppm_df.loc[i, "osm_method"],
            "Technology",
        ] = osm_ppm_df.loc[i, "ppm_technology"]

    # originates from osm::"tags.generator:source"
    add_ppls.loc[add_ppls["Fueltype"] == "Nuclear", "Technology"] = "Steam Turbine"

    # PMM contains data on NG, batteries and hydro storages
    # trying to catch some of them...
    # originates from osm::"tags.generator:source"
    add_ppls.loc[add_ppls["Fueltype"] == "battery", "Set"] = "Store"
    # originates from osm::tags.generator:type
    add_ppls.loc[add_ppls["Technology"] == "battery storage", "Set"] = "Store"

    add_ppls = add_ppls.replace(dict(Fueltype={"battery": "Other"})).drop(
        columns=["tags.generator:method", "geometry", "Area", "id"],
        errors="ignore",
    )

    to_csv_nafix(add_ppls, filepath_ppl_pm, index=False)

    return add_ppls


def add_custom_powerplants(
    ppl: pd.DataFrame,
    custom_powerplants_files: list[str],
    custom_powerplants_config: dict,
) -> pd.DataFrame:
    """
    Apply the configured method to each custom powerplant file.

    Countries affected by each method are determined from the ``Country``
    column of the corresponding file.

    Parameters
    ----------
    ppl : pd.DataFrame
        Powerplantmatching dataframe.
    custom_powerplants_files : list[str]
        Paths to custom powerplants CSV files.
    custom_powerplants_config : dict
        Configuration containing ``filepaths`` and ``method``. ``method`` may
        be a scalar applied to all files or a list matching ``filepaths``.

    Returns
    -------
    pd.DataFrame
        Powerplants dataframe including the configured custom powerplants.
    """
    methods = custom_powerplants_config.get("method", False)

    if not isinstance(methods, list):
        methods = [methods] * len(custom_powerplants_files)

    if len(methods) != len(custom_powerplants_files):
        raise ValueError(
            "The number of custom powerplant methods must match the number of filepaths."
        )

    allowed_methods = {False, "merge", "replace"}
    invalid_methods = [method for method in methods if method not in allowed_methods]

    if invalid_methods:
        raise ValueError(
            "Custom powerplant methods must be false, 'merge', or 'replace'; "
            f"found {invalid_methods}."
        )

    result = ppl.copy()

    for filepath, method in zip(custom_powerplants_files, methods):
        if not method:
            continue

        custom_ppls = read_csv_nafix(
            filepath,
            index_col=0,
            dtype={"bus": "str"},
        )

        if "Country" not in custom_ppls.columns:
            raise ValueError(
                f"Custom powerplant file '{filepath}' has no 'Country' column."
            )

        countries = custom_ppls["Country"].dropna().unique()

        if len(countries) == 0:
            raise ValueError(
                f"Custom powerplant file '{filepath}' contains no countries."
            )

        if method == "replace":
            result = result.loc[~result["Country"].isin(countries)]

        result = pd.concat(
            [result, custom_ppls],
            ignore_index=True,
            sort=False,
        )

    return result


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
        df.loc[fueltype, "Technology"] = df.loc[fueltype, "Technology"].replace(
            {t: "CCGT" for t in unknown_techs}
        )
    df["Fueltype"] = np.where(fueltype, df["Technology"], df["Fueltype"])
    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_powerplants")

    configure_logging(snakemake)

    with open(snakemake.input.pm_config, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    filepath_osm_ppl = snakemake.input.osm_powerplants
    filepath_osm2pm_ppl = snakemake.output.powerplants_osm2pm

    n = pypsa.Network(snakemake.input.base_network)
    countries_codes = n.buses.country.unique()
    countries_names = list(map(two_digits_2_name_country, countries_codes))

    config["target_countries"] = countries_names

    if (
        "EXTERNAL_DATABASE"
        in config["matching_sources"] + config["fully_included_sources"]
    ):
        if "EXTERNAL_DATABASE" not in config:
            logger.error(
                "Missing configuration EXTERNAL_DATABASE in powerplantmatching config yaml\n\t"
                "Please check file configs/powerplantmatching_config.yaml"
            )
        logger.info("Parsing OSM generator data to powerplantmatching format")
        config["EXTERNAL_DATABASE"]["fn"] = os.path.join(
            os.getcwd(), filepath_osm2pm_ppl
        )
    else:
        # create an empty file
        with open(filepath_osm2pm_ppl, "w"):
            pass

    # specify the main query for filtering powerplants
    ppl_query = snakemake.params.powerplants_filter
    if isinstance(ppl_query, str):
        config["main_query"] = ppl_query
    else:
        config["main_query"] = ""

    custom_powerplants = snakemake.params.custom_powerplants

    ppl = (
        pm.powerplants(from_url=False, update=True, config_update=config)
        .powerplant.fill_missing_decommissioning_years()
        .query("Country in @countries_names")
        .powerplant.convert_country_to_alpha2()
        .pipe(replace_natural_gas_technology)
    )

    ppl = add_custom_powerplants(
        ppl=ppl,
        custom_powerplants_files=list(snakemake.input.custom_powerplants),
        custom_powerplants_config=custom_powerplants,
    )

    if isinstance(ppl_query, str) and ppl_query:
        ppl = ppl.query(ppl_query)

    # define unique index
    ppl = ppl.reset_index(drop=True)

    cntries_without_ppl = [c for c in countries_codes if c not in ppl.Country.unique()]

    for c in countries_codes:
        substation_i = n.buses.query("substation_lv and country == @c").index
        kdtree = KDTree(n.buses.loc[substation_i, ["x", "y"]].values)
        ppl_i = ppl.query("Country == @c").index

        tree_i = kdtree.query(ppl.loc[ppl_i, ["lon", "lat"]].values)[1]
        ppl.loc[ppl_i, "bus"] = substation_i.append(pd.Index([np.nan]))[tree_i]

    if cntries_without_ppl:
        logger.warning(f"No powerplants known in: {', '.join(cntries_without_ppl)}")

    bus_null_b = ppl["bus"].isnull()
    if bus_null_b.any():
        logger.warning(f"Couldn't find close bus for {bus_null_b.sum()} powerplants")

    if snakemake.params.alternative_clustering:
        gadm_layer_id = snakemake.params.gadm_layer_id
        country_list = snakemake.params.countries
        geo_crs = snakemake.params.geo_crs

        ppl = locate_bus(
            ppl.rename(columns={"lon": "x", "lat": "y", "Country": "country"}),
            country_list,
            gadm_layer_id,
            snakemake.input.gadm_shapes,
            snakemake.params.alternative_clustering,
            col_out="region_id",
        ).rename(columns={"x": "lon", "y": "lat", "country": "Country"})

    ppl.to_csv(snakemake.output.powerplants)
