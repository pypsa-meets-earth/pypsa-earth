# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Retrieves conventional powerplant capacities and locations from `powerplantmatching <https://github.com/FRESNA/powerplantmatching>`_, assigns these to buses and creates a ``.csv`` file. It is possible to merge the powerplant database with, or replace it using, one or more custom powerplant files.

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
- Files listed under ``custom_powerplants.filepaths``: custom powerplants in the same format as `powerplantmatching <https://github.com/FRESNA/powerplantmatching>`_ provides or as the OSM extractor generates.

Outputs
-------

- ``resource/powerplants.csv``: A list of conventional power plants (i.e. neither wind nor solar) with fields for name, fuel type, technology, country, capacity in MW, duration, commissioning year, retrofit year, latitude, longitude, and dam information as documented in the `powerplantmatching README <https://github.com/FRESNA/powerplantmatching/blob/master/README.md>`_; additionally it includes information on the closest substation/bus in ``networks/base.nc``.

    .. image:: /img/powerplantmatching.png
        :width: 30 %

    **Source:** `powerplantmatching on GitHub <https://github.com/FRESNA/powerplantmatching>`_

Description
-----------

The configuration option ``electricity: powerplants_filter`` specifies a `pandas.query <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.query.html>`_ command applied to the original powerplantmatching database. Country conditions must therefore use the full country names found in the original database, such as ``United States``. Country values are converted to ISO alpha-2 codes after the custom powerplant files have been applied.

The ``electricity: custom_powerplants`` section specifies the custom files and how they are applied. ``method`` accepts ``false``, ``merge``, or ``replace`` and is applied once to the complete set of configured files. ``merge`` appends all configured custom powerplants to the filtered powerplantmatching dataset, while ``replace`` discards the complete powerplantmatching dataset and uses all configured custom files instead.

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

3. Replacing the complete powerplantmatching dataset:

    .. code:: yaml

        custom_powerplants:
          filepaths:
          - data/custom_powerplants.csv
          method: replace

4. Combining multiple custom files:

   Multiple custom powerplant files can be provided. A single ``method`` is
   applied to the complete set of configured files.

   .. code:: yaml

       custom_powerplants:
         filepaths:
         - data/custom_powerplants_US.csv
         - data/custom_powerplants_CA.csv
         method: merge

   The available methods are:

   - ``false``: ignore the custom files and use only powerplantmatching data.
   - ``merge``: add all configured custom files to the filtered
     powerplantmatching dataset.
   - ``replace``: discard the complete powerplantmatching dataset and use only
     the configured custom files.

5. Obtaining different outcomes for different countries:

   Country-specific replacement can be achieved by excluding the corresponding country from ``powerplants_filter`` and then merging the custom files.

   .. code:: yaml

       powerplants_filter: (DateOut >= 2022 or DateOut != DateOut) and (DateIn <= 2023 or DateIn != DateIn) and Country != 'United States'

       custom_powerplants:
         filepaths:
         - data/custom_powerplants_US.csv
         - data/custom_powerplants_CA.csv
         method: merge

   In this example:

   - US plants are taken only from ``custom_powerplants_US.csv`` because US
     plants are excluded from powerplantmatching.
   - CA plants from ``custom_powerplants_CA.csv`` are added to the existing
     powerplantmatching data.
   - Countries without a custom file, such as MX, use only powerplantmatching
     data.

Custom powerplant file format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Custom powerplant CSV files should follow the powerplantmatching format, with some additional considerations.

Required columns:

``id``, ``Name``, ``Fueltype``, ``Technology``, ``Set``, ``Country``,
``Capacity``, ``Efficiency``, ``DateIn``, ``DateRetrofit``, ``DateOut``,
``lat``, ``lon``, ``Duration``, ``Volume_Mm3``, ``DamHeight_m``,
``StorageCapacity_MWh``, ``EIC``, and ``projectID``.

Considerations for column values:

- ``Fueltype`` should use the corresponding powerplantmatching fuel category, such as ``Natural Gas``.
- Natural-gas plants should specify ``OCGT`` or ``CCGT`` in the ``Technology`` column.
- Hydro plants with ``Technology`` set to ``Reservoir`` are represented as storage units. Use ``ror`` for hydro plants that should be represented as generators.
- ``Country`` values are converted to ISO 3166-1 alpha-2 codes after the custom files have been applied. Both full country names and alpha-2 codes are therefore accepted in custom files.

OSM mapping assumptions
~~~~~~~~~~~~~~~~~~~~~~~

The following assumptions were made when mapping custom OSM-extracted power plants to the powerplantmatching format:

1. The benchmark powerplantmatching values were taken as follows:

   - ``Fueltype``: ``Hydro``, ``Hard Coal``, ``Natural Gas``, ``Lignite``,
     ``Nuclear``, ``Oil``, ``Bioenergy``, ``Wind``, ``Geothermal``, ``Solar``,
     ``Waste``, and ``Other``.
   - ``Technology``: ``Reservoir``, ``Pumped Storage``, ``Run-Of-River``,
     ``Steam Turbine``, ``CCGT``, ``OCGT``, ``PV``, ``CCGT, Thermal``,
     ``Offshore``, and ``Storage Technologies``.
   - ``Set``: ``Store``, ``PP``, and ``CHP``.

2. OSM-extracted features were mapped to powerplantmatching values using the following rules:

   - ``coal`` -> ``Hard Coal``
   - ``wind_turbine`` -> ``Onshore``
   - ``horizontal_axis`` -> ``Onshore``
   - ``vertical_axis`` -> ``Offshore``
   - ``nuclear`` -> ``Steam Turbine``

3. All hydro objects extracted from OSM were interpreted as generation technologies, although ``Run-Of-River``, ``Pumped Storage``, and ``Reservoir`` may also belong to ``Storage Technologies`` in powerplantmatching.

4. The OSM extraction was assumed to ignore non-generation features such as CHP plants and natural-gas storage, unlike powerplantmatching.
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
    method: bool | str,
) -> pd.DataFrame:
    """
    Merge or replace powerplantmatching data with custom powerplant files.

    Parameters
    ----------
    ppl : pd.DataFrame
        Powerplantmatching dataframe.
    custom_powerplants_files : list[str]
        Paths to custom powerplant CSV files.
    method : bool or str
        Method applied to all custom files together. Accepted values are
        ``False``, ``"merge"``, and ``"replace"``.

    Returns
    -------
    pd.DataFrame
        Powerplant dataframe after applying the custom files.
    """
    allowed_methods = {False, "merge", "replace"}
    if method not in allowed_methods:
        raise ValueError(
            "electricity.custom_powerplants.method must be false, "
            f"'merge', or 'replace'; found {method!r}."
        )

    if method is False:
        return ppl

    if not custom_powerplants_files:
        raise ValueError(
            "At least one filepath must be configured when "
            "electricity.custom_powerplants.method is 'merge' or 'replace'."
        )

    custom_ppls = pd.concat(
        [
            read_csv_nafix(
                filepath,
                index_col=0,
                dtype={"bus": "str"},
            )
            for filepath in custom_powerplants_files
        ],
        ignore_index=True,
        sort=False,
    )

    if method == "replace":
        return custom_ppls

    return pd.concat(
        [ppl, custom_ppls],
        ignore_index=True,
        sort=False,
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
    custom_method = custom_powerplants["method"]

    if custom_method == "replace":
        ppl = pd.DataFrame()
    else:
        ppl = (
            pm.powerplants(
                from_url=False,
                update=True,
                config_update=config,
            )
            .powerplant.fill_missing_decommissioning_years()
            .query("Country in @countries_names")
        )

    ppl = add_custom_powerplants(
        ppl=ppl,
        custom_powerplants_files=list(snakemake.input.custom_powerplants),
        method=custom_method,
    )

    ppl = ppl.powerplant.convert_country_to_alpha2().pipe(
        replace_natural_gas_technology
    )

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
