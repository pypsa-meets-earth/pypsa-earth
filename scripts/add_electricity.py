# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
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

    renewable:
        hydro:
            carriers:
            hydro_max_hours:
            hydro_capital_cost:

    lines:
        length_factor:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at :ref:`costs_cf`,
    :ref:`electricity_cf`, :ref:`load_options_cf`, :ref:`renewable_cf`, :ref:`lines_cf`

Inputs
------

- ``resources/costs.csv``: The database of cost assumptions for all included technologies for specific years from various sources; e.g. discount rate, lifetime, investment (CAPEX), fixed operation and maintenance (FOM), variable operation and maintenance (VOM), fuel costs, efficiency, carbon-dioxide intensity.
- ``data/bundle/hydro_capacities.csv``: Hydropower plant store/discharge power capacities, energy storage capacity, and average hourly inflow by country.  Not currently used!

    .. image:: /img/hydrocapacities.png

- ``data/geth2015_hydro_capacities.csv``: alternative to capacities above; not currently used!
- ``resources/demand_profiles.csv``: a csv file containing the demand profile associated with buses
- ``resources/shapes/gadm_shapes.geojson``: confer :ref:`shapes`
- ``resources/powerplants.csv``: confer :ref:`powerplants`
- ``resources/profile_{}.nc``: all technologies in ``config["renewables"].keys()``, confer :ref:`renewableprofiles`
- ``networks/base.nc``: confer :ref:`base`

Outputs
-------

- ``networks/elec.nc``:

    .. image:: /img/elec.png
            :width: 75 %
            :align: center

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
from typing import Dict, List

import geopandas as gpd
import numpy as np
import pandas as pd
import powerplantmatching as pm
import pypsa
import xarray as xr
from _helpers import configure_logging, getContinent, update_p_nom_max
from config_osm_data import world_iso
from powerplantmatching.export import map_country_bus
from shapely.validation import make_valid
from vresutils import transfer as vtransfer

idx = pd.IndexSlice

logger = logging.getLogger(__name__)

carrier_pypsa_mapping = {
    "ocgt": "OCGT",
    "ccgt": "CCGT",
    "bioenergy": "biomass",
    "ccgt, thermal": "CCGT",
    "hard coal": "coal",
}
# rename because technology data & pypsa earth costs.csv use different names
# TODO: rename the technologies in hosted tutorial data to match technology data

hydrogen_pypsa_mapping = {
    "hydrogen storage": "hydrogen storage tank",
    "hydrogen storage tank": "hydrogen storage tank",
    "hydrogen storage tank type 1": "hydrogen storage tank",
    "hydrogen underground storage": "hydrogen storage underground",
}


def normed(s):
    return s / s.sum()


def calculate_annuity(n, r):
    """
    Calculate the annuity factor for an asset with lifetime n years and
    discount rate of r, e.g. annuity(20, 0.05) * 20 = 1.6
    """
    if isinstance(r, pd.Series):
        return pd.Series(1 / n, index=r.index).where(
            r == 0, r / (1.0 - 1.0 / (1.0 + r) ** n)
        )
    elif r > 0:
        return r / (1.0 - 1.0 / (1.0 + r) ** n)
    else:
        return 1 / n


def _add_missing_carriers_from_costs(n, costs, carriers, countries):
    for country in countries:
        country_mask = costs.index.get_level_values("country") == country
        missing_carriers = pd.Index(carriers).difference(n.carriers.index)
        if missing_carriers.empty:
            continue
        costs_filtered = costs.loc[country]
        emissions_cols = (
            costs_filtered.columns.to_series()
            .loc[lambda s: s.str.endswith("_emissions")]
            .values
        )
        suptechs = missing_carriers.str.split("-").str[0]
        emissions = costs_filtered.loc[suptechs, emissions_cols].fillna(0.0)

        emissions.index = missing_carriers
        emissions.pipe(
            lambda x: x.assign(country=country)
            .set_index("country", append=True)
            .swaplevel(0, 1)
        )
        n.import_components_from_dataframe(emissions, "Carrier")


def load_costs(
    tech_costs_path: str,
    config_costs: Dict,
    elec_config: Dict,
    countries: List[str],
    country_mapping: Dict[str, str],
    Nyears: int = 1,
) -> pd.DataFrame:
    """
    Set all asset costs and formats them according to preferences.

    Parameters
    ----------
    tech_costs_path : str
        Path to the technology costs file.
    config_costs : Dict
        Dictionary of costs section defined in the config.yaml.
    elec_config : Dict
        Dictionary of electricity section defined in the config.yaml.
    countries : List[str]
        List of country ISO2 codes based on the countries defined in the config.
    country_mapping : Dict[str,str]
       Dictionary of country ISO2 keys organised by continent.
    Nyears : int, optional
        _description_, by default 1

    Returns
    -------
    pd.DataFrame
        Processed cost data.

    """
    max_hours = elec_config["max_hours"]

    costs = (
        pd.read_csv(tech_costs_path, index_col=["country", "technology", "parameter"])
        .sort_index()
        .pipe(correct_units, config_costs)
        .pipe(fill_nas_with_config_costs, config_costs)
        .pipe(annualise_capital_costs, Nyears)
        .rename(columns={"CO2 intensity": "co2_emissions"})
        .pipe(transform_country_costs, countries, country_mapping)
        .pipe(update_gas_costs, countries)
        .pipe(calculate_solar_capital_costs, config_costs, countries)
        .pipe(apply_countrywise_battery_costs, max_hours, countries)
        .pipe(apply_overwrites_from_config, countries, config_costs)
    )

    return costs


def load_global_costs(
    tech_costs_path: str,
    Nyears: float,
    config_costs: Dict,
    countries: List[str] = ["global"],
) -> pd.DataFrame:
    """Specific function to only retrieve global costs. Used for simplify network cost calculations.

    Parameters
    ----------
    tech_costs_path : str
        Path to the technology costs file.
    Nyears : float
        _description_
    config_costs : Dict
        Dictionary of costs section defined in the config.yaml.
    countries : List[str], optional
        list with global., by default ["global"]

    Returns
    -------
    pd.DataFrame
       cost database specifically for calculating offshore wind costs in simplify network.
    """
    costs = (
        pd.read_csv(tech_costs_path, index_col=["country", "technology", "parameter"])
        .sort_index()
        .pipe(correct_units, config_costs)
        .pipe(fill_nas_with_config_costs, config_costs)
        .pipe(annualise_capital_costs, Nyears)
        .rename(columns={"CO2 intensity": "co2_emissions"})
        .pipe(update_gas_costs, countries)
        .pipe(calculate_solar_capital_costs, config_costs, countries)
        .pipe(apply_overwrites_from_config, countries, config_costs)
    )

    return costs


def rename_technologies(costs: pd.DataFrame, mapping: Dict[str, str]) -> pd.DataFrame:
    costs = costs.rename(mapping)

    return costs


def apply_overwrites_from_config(
    costs: pd.DataFrame, countries: List[str], config_costs: Dict
) -> pd.DataFrame:
    """Applies overwrites from the config file for marginal and capital costs if defined.

    Parameters
    ----------
    costs : pd.DataFrame
         Technology-costs dataframe.
    countries : List[str]
       List of country ISO2 codes based on the countries defined in the config.
    config_costs : Dict
        Dictionary of costs section defined in the config.yaml.

    Returns
    -------
    pd.DataFrame
       Technology costs dataframe with overwrites applied.
    """
    for country in countries:
        country_mask = costs.index.get_level_values("country") == country
        # TODO  check that lines below are functioning as desired
        for attr in ("marginal_cost", "capital_cost"):
            overwrites = config_costs.get(attr)
            if overwrites is not None:
                overwrites = pd.Series(overwrites)
                costs.loc[country_mask, attr] = overwrites.get(
                    country, costs.loc[country_mask, attr]
                )
    return costs


def apply_countrywise_battery_costs(
    costs: pd.DataFrame, max_hours: Dict, countries: List[str]
) -> pd.DataFrame:
    """Applies battery costs by country.

    Parameters
    ----------
    costs : pd.DataFrame
        Technology-costs dataframe.
    max_hours : Dict
        _description_
    countries : List[str]
        List of country ISO2 codes based on the countries defined in the config.

    Returns
    -------
    pd.DataFrame
       Costs dataframe with country specific battery costs applied.
    """
    for country in countries:
        battery_storage_cost = costs.loc[(country, "battery storage")]
        battery_inverter_cost = costs.loc[(country, "battery inverter")]
        battery_cost = costs_for_storage(
            battery_storage_cost, battery_inverter_cost, max_hours=max_hours["battery"]
        )
        costs.loc[(country, "battery")] = battery_cost
    return costs


def correct_units(costs: pd.DataFrame, config_costs: Dict) -> pd.DataFrame:
    """Corrects units to MW and EUR

    Parameters
    ----------
    costs : pd.DataFrame
        Technology-costs dataframe.
    config_costs : Dict
        Dictionary of costs section defined in the config.yaml.

    Returns
    -------
    pd.DataFrame
        _description_
    """
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3
    costs.unit = costs.unit.str.replace("/kW", "/MW")
    costs.loc[costs.unit.str.contains("USD"), "value"] *= config_costs[
        "USD2013_to_EUR2013"
    ]  # TODO #777 generalise currency conversion

    return costs


def fill_nas_with_config_costs(costs: pd.DataFrame, config_costs: Dict) -> pd.DataFrame:
    """Fills any NA values in costs with the fill_values defined in the config.

    Parameters
    ----------
    costs : pd.DataFrame
       Technology-costs dataframe.
    config_costs : Dict
        Dictionary of costs section defined in the config.yaml.

    Returns
    -------
    pd.DataFrame
        costs dataframe updated with config fill values.
    """
    costs = costs.value.unstack().fillna(config_costs["fill_values"])
    return costs


def annualise_capital_costs(costs: pd.DataFrame, Nyears: float) -> pd.DataFrame:
    """Annualises capital costs.

    Parameters
    ----------
    costs : pd.DataFrame
        Technology-costs dataframe.
    Nyears : float
        _description_

    Returns
    -------
    pd.DataFrame
        costs dataframe with annualised capital costs.
    """

    costs["capital_cost"] = (
        (
            calculate_annuity(costs["lifetime"], costs["discount rate"])
            + costs["FOM"] / 100.0
        )
        * costs["investment"]
        * Nyears
    )
    return costs


def update_gas_costs(costs: pd.DataFrame, countries: List[str]) -> pd.DataFrame:
    """
    Updates gas costs and emissions based on the gas technology and the country of origin.

    Parameters
    ----------
    costs : pd.DataFrame
        Technology-costs dataframe.
    countries : List[str]
        List of countries defined in the config.

    Returns
    -------
    pd.DataFrame
         Updated costs dataframe with solar capital costs.
    """
    for country in countries:
        country_mask = costs.index.get_level_values("country") == country
        mask_ocgt = (costs.index.get_level_values("country") == country) & (
            costs.index.get_level_values("technology") == "OCGT"
        )
        mask_ccgt = (costs.index.get_level_values("country") == country) & (
            costs.index.get_level_values("technology") == "CCGT"
        )
        mask_gas = (costs.index.get_level_values("country") == country) & (
            costs.index.get_level_values("technology") == "gas"
        )
        costs.loc[mask_ocgt, "fuel"] = costs.loc[mask_gas, "fuel"].values
        costs.loc[mask_ocgt, "co2_emissions"] = costs.loc[
            mask_gas, "co2_emissions"
        ].values
        costs.loc[mask_ccgt, "fuel"] = costs.loc[mask_gas, "fuel"].values
        costs.loc[mask_ccgt, "co2_emissions"] = costs.loc[
            mask_gas, "co2_emissions"
        ].values
        costs.loc[country_mask, "marginal_cost"] = (
            costs.loc[country_mask, "VOM"]
            + costs.loc[country_mask, "fuel"] / costs.loc[country_mask, "efficiency"]
        )

    return costs


def calculate_solar_capital_costs(
    costs: pd.DataFrame, config_costs: Dict, countries: List[str]
) -> pd.DataFrame:
    """
    Calculates capital costs for solar-rooftop and solar-utility based on the country of origin and the share of each technology

    Parameters
    ----------
    costs : pd.DataFrame
    Technology-costs dataframe.
    config_costs : Dict
        Dictionary of costs section defined in the config.yaml.
    countries : List[str]
        List of countries defined in the config.

    Returns
    -------
    pd.DataFrame
        Updated costs dataframe with solar capital costs.
    """
    for country in countries:
        rooftop_share = config_costs.get(country, {}).get("rooftop_share", 0.0)
        country_mask = costs.index.get_level_values("country") == country
        solar_rooftop_cost = costs.loc[(country, "solar-rooftop"), "capital_cost"]
        solar_utility_cost = costs.loc[(country, "solar-utility"), "capital_cost"]
        utility_share = 1 - rooftop_share
        costs.loc[
            country_mask & (costs.index.get_level_values("technology") == "solar"),
            "capital_cost",
        ] = (
            rooftop_share * solar_rooftop_cost + utility_share * solar_utility_cost
        )
    return costs


def costs_for_storage(store, link1, link2=None, max_hours=1.0):
    # TODO reduce verbosity and duplication
    capital_cost = link1["capital_cost"] + max_hours * store["capital_cost"]
    if link2 is not None:
        capital_cost += link2["capital_cost"]
    return pd.Series(
        dict(capital_cost=capital_cost, marginal_cost=0.0, co2_emissions=0.0)
    )


def load_powerplants(file_path: str) -> pd.DataFrame:
    """
    Load power plant data from a CSV file and perform data transformations.

    Parameters
    ----------
    file_path : str
        Path to the CSV file containing power plant data.

    carrier_mapping : Dict[str, str]
        A dictionary mapping carrier values to their replacements based on technology data.

    Returns
    -------
    pd.DataFrame
        Processed power plant data.

    Raises
    ------
    FileNotFoundError
        If the specified file_path does not exist.
    """

    try:
        powerplants = pd.read_csv(file_path, index_col=0, dtype={"bus": "str"})
    except FileNotFoundError:
        raise FileNotFoundError("File not found: {}".format(file_path))

    powerplants = (
        powerplants.powerplant.to_pypsa_names()
        .powerplant.convert_country_to_alpha2()
        .rename(columns=str.lower)
        .drop(columns=["efficiency"])
    )
    return powerplants


def map_carriers(
    powerplants: pd.DataFrame, carrier_mapping: Dict[str, str]
) -> pd.DataFrame:
    powerplants = powerplants.replace({"carrier": carrier_mapping})

    return powerplants


def attach_load(n: pypsa.Network, demand_profiles_path: str) -> None:
    """
    Add load profiles to network buses.

    Parameters
    ----------
    n: pypsa.Network
        The PyPSA network object.

    demand_profiles: str
        Path to csv file of elecric demand time series, e.g. "resources/demand_profiles.csv"
        Demand profile has snapshots as rows and bus names as columns.


    """
    demand_df = pd.read_csv(demand_profiles_path, index_col=0, parse_dates=True)

    n.madd("Load", demand_df.columns, bus=demand_df.columns, p_set=demand_df)


def update_transmission_costs(
    n, costs, countries: List[str], length_factor=1.0, simple_hvdc_costs=False
):
    lines_mask = n.lines["country"].isin(countries)
    lines_filtered = n.lines.loc[lines_mask]

    hvac_overhead_costs = costs.loc[
        (countries, "HVAC overhead"), "capital_cost"
    ].to_dict()
    hvdc_overhead_costs = costs.loc[
        (countries, "HVDC overhead"), "capital_cost"
    ].to_dict()

    for country in countries:
        country_mask = lines_filtered["country"] == country
        hvac_overhead_cost = hvac_overhead_costs.get(country, 0.0)

        n.lines.loc[country_mask, "capital_cost"] = (
            lines_filtered.loc[country_mask, "length"].values
            * length_factor
            * hvac_overhead_cost
        )
    if n.links.empty:
        return

    dc_mask = (n.links["carrier"] == "DC") & n.links["country"].isin(countries)
    links_dc_filtered = n.links.loc[dc_mask]
    # If there are no "DC" links, then the 'underwater_fraction' column
    # may be missing. Therefore we have to return here.
    # TODO: Require fix
    if links_dc_filtered.empty:
        return

    if simple_hvdc_costs:
        hdvc_overhead_costs = costs.loc[
            (countries, "HVDC overhead"), "capital_cost"
        ].to_dict()

        n.links.loc[dc_mask, "capital_cost"] = (
            links_dc_filtered["length"].values * length_factor * hdvc_overhead_costs
        )
    else:
        underwater_fraction = links_dc_filtered["underwater_fraction"].values
        hvdc_submarine_costs = costs.loc[
            (
                countries,
                "HVDC submarine",
            ),
            "capital_cost",
        ].to_dict()
        hvdc_inverter_pair_costs = costs.loc[
            (countries, "HVDC inverter pair"), "capital_cost"
        ].to_dict()

        n.links.loc[dc_mask, "capital_cost"] = links_dc_filtered[
            "length"
        ].values * length_factor * (
            (1.0 - underwater_fraction) * hvdc_overhead_costs.get(country, 0.0)
            + underwater_fraction * hvdc_submarine_costs.get(country, 0.0)
        ) + hvdc_inverter_pair_costs.get(
            country, 0.0
        )


def attach_wind_and_solar(
    n,
    costs,
    ppl,
    input_profiles,
    countries,
    technologies,
    extendable_carriers,
    line_length_factor=1,
):
    # TODO: rename tech -> carrier, technologies -> carriers
    _add_missing_carriers_from_costs(n, costs, technologies, countries)

    df = ppl.rename(columns={"country": "Country"})

    for tech in technologies:
        costs_filtered = costs.loc[
            costs.index.get_level_values("technology") == tech
        ].droplevel("technology")

        if tech == "hydro":
            continue

        if tech == "offwind-ac":
            # add all offwind wind power plants by default as offwind-ac
            df.carrier.mask(df.technology == "Offshore", "offwind-ac", inplace=True)

        df.carrier.mask(df.technology == "Onshore", "onwind", inplace=True)

        with xr.open_dataset(getattr(snakemake.input, "profile_" + tech)) as ds:
            if ds.indexes["bus"].empty:
                continue
            geo_costs = n.buses.loc[ds.indexes["bus"]].country
            costs_by_bus = pd.DataFrame(geo_costs).join(
                costs_filtered, on="country", rsuffix="_r"
            )
            suptech = tech.split("-", 2)[0]
            if suptech == "offwind":
                continue
                # TODO: Uncomment out and debug.
                # underwater_fraction = ds["underwater_fraction"].to_pandas()
                # connection_cost = (
                #     snakemake.config["lines"]["length_factor"] *
                #     ds["average_distance"].to_pandas() *
                #     (underwater_fraction *
                #      costs.at[tech + "-connection-submarine", "capital_cost"] +
                #      (1.0 - underwater_fraction) *
                #      costs.at[tech + "-connection-underground", "capital_cost"]
                #      ))
                # capital_cost = (costs.at["offwind", "capital_cost"] +
                #                 costs.at[tech + "-station", "capital_cost"] +
                #                 connection_cost)
                # logger.info(
                #     "Added connection cost of {:0.0f}-{:0.0f} Eur/MW/a to {}".
                #     format(connection_cost.min(), connection_cost.max(), tech))
            if not df.query("carrier == @tech").empty:
                buses = n.buses.loc[ds.indexes["bus"]]
                caps = map_country_bus(df.query("carrier == @tech"), buses)
                caps = caps.groupby(["bus"]).p_nom.sum()
                caps = pd.Series(data=caps, index=ds.indexes["bus"]).fillna(0)
            else:
                caps = pd.Series(index=ds.indexes["bus"]).fillna(0)

            n.madd(
                "Generator",
                ds.indexes["bus"],
                " " + tech,
                bus=ds.indexes["bus"],
                carrier=tech,
                p_nom=caps,
                p_nom_extendable=tech in extendable_carriers["Generator"],
                p_nom_min=caps,
                p_nom_max=ds["p_nom_max"].to_pandas(),
                p_max_pu=ds["profile"].transpose("time", "bus").to_pandas(),
                weight=ds["weight"].to_pandas(),
                marginal_cost=costs_by_bus["marginal_cost"].values,
                capital_cost=costs_by_bus["capital_cost"].values,
                efficiency=costs_by_bus["efficiency"].values,
            )


def transform_country_costs(
    costs: pd.DataFrame, countries: List[str], country_mapping: Dict[str, str]
) -> pd.DataFrame:
    """Transforms global or regional index strings into relevant country ones.

    Parameters
    ----------
    costs : pd.DataFrame
        Pandas dataframe of costs.
    countries : List[str]
        List of country ISO2 codes based on the countries defined in the config.
    country_mapping : Dict[str,str]
        Dictionary of country ISO2 keys organised by continent.
    """

    cost_index = costs.index.get_level_values("country")
    transformed_costs_dfs = []
    for country in countries:
        continent = find_parent_key(country_mapping, country)
        if country in cost_index:
            country_costs = costs.loc[country].assign(country=country)
        elif continent in cost_index:
            country_costs = costs.loc[continent].assign(country=country)
        else:
            country_costs = costs.loc["global"].assign(country=country)

        transformed_costs_dfs.append(country_costs)

    transformed_df = (
        pd.concat(transformed_costs_dfs)
        .set_index("country", append=True)
        .swaplevel(0, 1)
    )
    transformed_df = transformed_df.loc[
        transformed_df.index.isin(countries, level="country")
    ]

    return transformed_df


def find_parent_key(dictionary: Dict[str, str], target_key: str):
    """Looks for country key inside continent.

    Parameters
    ----------
    dictionary : Dict[str, str]
        dictionary of continents and countries
    target_key : str
        country iso2

    Returns
    -------
        continent key
    """
    for parent_key, child_dict in dictionary.items():
        if target_key in child_dict:
            return parent_key
    return None  # Return None if the target key is not found


def attach_conventional_generators(
    n,
    costs,
    ppl,
    conventional_carriers,
    extendable_carriers,
    renewable_carriers,
    conventional_config,
    conventional_inputs,
    countries,
):
    carriers = set(conventional_carriers) | set(extendable_carriers["Generator"]) - set(
        renewable_carriers
    )
    _add_missing_carriers_from_costs(n, costs, carriers, countries)

    ppl = (
        ppl.query("carrier in @carriers")
        .join(costs, on=["country", "carrier"], rsuffix="_r")
        .rename(index=lambda s: "C" + str(s))
    )
    ppl["efficiency"] = ppl.efficiency.fillna(ppl.efficiency)

    logger.info(
        "Adding {} generators with capacities [GW] \n{}".format(
            len(ppl), ppl.groupby("carrier").p_nom.sum().div(1e3).round(2)
        )
    )

    n.madd(
        "Generator",
        ppl.index,
        carrier=ppl.carrier,
        bus=ppl.bus,
        p_nom_min=ppl.p_nom.where(ppl.carrier.isin(conventional_carriers), 0),
        p_nom=ppl.p_nom.where(ppl.carrier.isin(conventional_carriers), 0),
        p_nom_extendable=ppl.carrier.isin(extendable_carriers["Generator"]),
        efficiency=ppl.efficiency,
        marginal_cost=ppl.marginal_cost,
        capital_cost=ppl.capital_cost,
        build_year=ppl.datein.fillna(0).astype(int),
        lifetime=(ppl.dateout - ppl.datein).fillna(np.inf),
    )

    for carrier in conventional_config:
        # Generators with technology affected
        idx = n.generators.query("carrier == @carrier").index

        for attr in list(set(conventional_config[carrier]) & set(n.generators)):
            values = conventional_config[carrier][attr]

            if f"conventional_{carrier}_{attr}" in conventional_inputs:
                # Values affecting generators of technology k country-specific
                # First map generator buses to countries; then map countries to p_max_pu
                values = pd.read_csv(values, index_col=0).iloc[:, 0]
                bus_values = n.buses.country.map(values)
                n.generators[attr].update(
                    n.generators.loc[idx].bus.map(bus_values).dropna()
                )
            else:
                # Single value affecting all generators of technology k indiscriminantely of country
                n.generators.loc[idx, attr] = values


def attach_hydro(n, costs, ppl, countries):
    if "hydro" not in snakemake.config["renewable"]:
        return
    c = snakemake.config["renewable"]["hydro"]
    carriers = c.get("carriers", ["ror", "PHS", "hydro"])

    _add_missing_carriers_from_costs(n, costs, carriers, countries)

    ppl = (
        ppl.query('carrier == "hydro"')
        .reset_index(drop=True)
        .rename(index=lambda s: str(s) + " hydro")
    ).join(costs, on=["country", "carrier"], rsuffix="_r")

    # TODO: remove this line to address nan when powerplantmatching is stable
    # Current fix, NaN technologies set to ROR
    ppl.loc[ppl.technology.isna(), "technology"] = "Run-Of-River"

    ror = ppl.query('technology == "Run-Of-River"')
    phs = ppl.query('technology == "Pumped Storage"')
    hydro = ppl.query('technology == "Reservoir"')

    if snakemake.config["cluster_options"]["alternative_clustering"]:
        bus_id = ppl["region_id"]
    else:
        bus_id = ppl["bus"]

    inflow_idx = ror.index.union(hydro.index)
    if not inflow_idx.empty:
        with xr.open_dataarray(snakemake.input.profile_hydro) as inflow:
            inflow_buses = bus_id[inflow_idx]
            missing_plants = pd.Index(inflow_buses.unique()).difference(
                inflow.indexes["plant"]
            )
            intersection_plants = pd.Index(
                inflow_buses[inflow_buses.isin(inflow.indexes["plant"])]
            )

            # if missing time series are found, notify the user and exclude missing hydro plants
            if not missing_plants.empty:
                # original total p_nom
                total_p_nom = ror.p_nom.sum() + hydro.p_nom.sum()
                idxs_to_keep = inflow_buses[
                    inflow_buses.isin(intersection_plants)
                ].index
                ror = ror.loc[ror.index.intersection(idxs_to_keep)]
                hydro = hydro.loc[hydro.index.intersection(idxs_to_keep)]
                # loss of p_nom
                loss_p_nom = ror.p_nom.sum() + hydro.p_nom.sum() - total_p_nom

                logger.warning(
                    f"'{snakemake.input.profile_hydro}' is missing inflow time-series for at least one bus: {', '.join(missing_plants)}."
                    f"Corresponding hydro plants are dropped, corresponding to a total loss of {loss_p_nom}MW out of {total_p_nom}MW."
                )

            if not intersection_plants.empty:
                inflow_t = (
                    inflow.sel(plant=intersection_plants)
                    .rename({"plant": "name"})
                    .assign_coords(name=inflow_idx)
                    .transpose("time", "name")
                    .to_pandas()
                )

    if "ror" in carriers and not ror.empty:
        n.madd(
            "Generator",
            ror.index,
            carrier="ror",
            bus=ror["bus"],
            p_nom=ror["p_nom"],
            efficiency=ror["efficiency"],
            capital_cost=ror["capital_cost"],
            weight=ror["p_nom"],
            p_max_pu=(
                inflow_t[ror.index]
                .divide(ror["p_nom"], axis=1)
                .where(lambda df: df <= 1.0, other=1.0)
            ),
        )

    if "PHS" in carriers and not phs.empty:
        # fill missing max hours to config value and
        # assume no natural inflow due to lack of data
        phs = phs.replace({"max_hours": {0: c["PHS_max_hours"]}})
        n.madd(
            "StorageUnit",
            phs.index,
            carrier="PHS",
            bus=phs["bus"],
            p_nom=phs["p_nom"],
            capital_cost=phs["capital_cost"],
            max_hours=phs["max_hours"],
            efficiency_store=np.sqrt(phs["efficiency"]),
            efficiency_dispatch=np.sqrt(phs["efficiency"]),
            cyclic_state_of_charge=True,
        )

    if "hydro" in carriers and not hydro.empty:
        hydro_max_hours = c.get("hydro_max_hours")
        hydro_stats = pd.read_csv(
            snakemake.input.hydro_capacities, comment="#", na_values=["-"], index_col=0
        )
        e_target = hydro_stats["E_store[TWh]"].clip(lower=0.2) * 1e6
        e_installed = hydro.eval("p_nom * max_hours").groupby(hydro.country).sum()
        e_missing = e_target - e_installed
        missing_mh_i = hydro.query("max_hours == 0").index

        if hydro_max_hours == "energy_capacity_totals_by_country":
            max_hours_country = (
                e_missing / hydro.loc[missing_mh_i].groupby("country").p_nom.sum()
            )

        elif hydro_max_hours == "estimate_by_large_installations":
            max_hours_country = (
                hydro_stats["E_store[TWh]"] * 1e3 / hydro_stats["p_nom_discharge[GW]"]
            )

        missing_countries = pd.Index(hydro["country"].unique()).difference(
            max_hours_country.dropna().index
        )
        if not missing_countries.empty:
            logger.warning(
                "Assuming max_hours=6 for hydro reservoirs in the countries: {}".format(
                    ", ".join(missing_countries)
                )
            )
        hydro_max_hours = hydro.max_hours.where(
            hydro.max_hours > 0, hydro.country.map(max_hours_country)
        ).fillna(6)

        n.madd(
            "StorageUnit",
            hydro.index,
            carrier="hydro",
            bus=hydro["bus"],
            p_nom=hydro["p_nom"],
            max_hours=hydro_max_hours,
            capital_cost=(
                costs.at["hydro", "capital_cost"]
                if c.get("hydro_capital_cost")
                else 0.0
            ),
            marginal_cost=hydro["marginal_cost"],
            p_max_pu=1.0,  # dispatch
            p_min_pu=0.0,  # store
            efficiency_dispatch=hydro["efficiency"],
            efficiency_store=0.0,
            cyclic_state_of_charge=True,
            inflow=inflow_t.loc[:, hydro.index],
        )


def attach_extendable_generators(n, costs, ppl):
    logger.warning("The function is deprecated with the next release")
    elec_opts = snakemake.config["electricity"]
    carriers = pd.Index(elec_opts["extendable_carriers"]["Generator"])

    _add_missing_carriers_from_costs(n, costs, carriers)

    for tech in carriers:
        if tech.startswith("OCGT"):
            ocgt = (
                ppl.query("carrier in ['OCGT', 'CCGT']")
                .groupby("bus", as_index=False)
                .first()
            )
            n.madd(
                "Generator",
                ocgt.index,
                suffix=" OCGT",
                bus=ocgt["bus"],
                carrier=tech,
                p_nom_extendable=True,
                p_nom=0.0,
                capital_cost=costs.at["OCGT", "capital_cost"],
                marginal_cost=costs.at["OCGT", "marginal_cost"],
                efficiency=costs.at["OCGT", "efficiency"],
            )

        elif tech.startswith("CCGT"):
            ccgt = (
                ppl.query("carrier in ['OCGT', 'CCGT']")
                .groupby("bus", as_index=False)
                .first()
            )
            n.madd(
                "Generator",
                ccgt.index,
                suffix=" CCGT",
                bus=ccgt["bus"],
                carrier=tech,
                p_nom_extendable=True,
                p_nom=0.0,
                capital_cost=costs.at["CCGT", "capital_cost"],
                marginal_cost=costs.at["CCGT", "marginal_cost"],
                efficiency=costs.at["CCGT", "efficiency"],
            )

        elif tech.startswith("nuclear"):
            nuclear = (
                ppl.query("carrier == 'nuclear'").groupby("bus", as_index=False).first()
            )
            n.madd(
                "Generator",
                nuclear.index,
                suffix=" nuclear",
                bus=nuclear["bus"],
                carrier=tech,
                p_nom_extendable=True,
                p_nom=0.0,
                capital_cost=costs.at["nuclear", "capital_cost"],
                marginal_cost=costs.at["nuclear", "marginal_cost"],
                efficiency=costs.at["nuclear", "efficiency"],
            )

        else:
            raise NotImplementedError(
                f"Adding extendable generators for carrier "
                "'{tech}' is not implemented, yet. "
                "Only OCGT, CCGT and nuclear are allowed at the moment."
            )


def estimate_renewable_capacities_irena(n, config):
    if not config["electricity"].get("estimate_renewable_capacities"):
        return

    stats = config["electricity"]["estimate_renewable_capacities"]["stats"]
    if not stats:
        return

    year = config["electricity"]["estimate_renewable_capacities"]["year"]
    tech_map = config["electricity"]["estimate_renewable_capacities"][
        "technology_mapping"
    ]
    tech_keys = list(tech_map.keys())
    countries = config["countries"]

    p_nom_max = config["electricity"]["estimate_renewable_capacities"]["p_nom_max"]
    p_nom_min = config["electricity"]["estimate_renewable_capacities"]["p_nom_min"]

    if len(countries) == 0:
        return
    if len(tech_map) == 0:
        return

    if stats == "irena":
        capacities = pm.data.IRENASTAT().powerplant.convert_country_to_alpha2()
    else:
        logger.info(
            f"Selected renewable capacity estimation statistics {stats} is not available, applying greenfield scenario instead"
        )
        return

    # Check if countries are in country list of stats
    missing = list(set(countries).difference(capacities.Country.unique()))
    if missing:
        logger.info(
            f"The countries {missing} are not provided in the stats and hence not scaled"
        )

    capacities = capacities.query(
        "Year == @year and Technology in @tech_keys and Country in @countries"
    )
    capacities = capacities.groupby(["Technology", "Country"]).Capacity.sum()

    logger.info(
        f"Heuristics applied to distribute renewable capacities [MW] "
        f"{capacities.groupby('Country').sum()}"
    )

    for ppm_technology, techs in tech_map.items():
        if ppm_technology not in capacities.index:
            logger.info(
                f"technology {ppm_technology} is not provided by {stats} and therefore not estimated"
            )
            continue

        tech_capacities = capacities.loc[ppm_technology].reindex(
            countries, fill_value=0.0
        )
        tech_i = n.generators.query("carrier in @techs").index
        n.generators.loc[tech_i, "p_nom"] = (
            (
                n.generators_t.p_max_pu[tech_i].mean()
                * n.generators.loc[tech_i, "p_nom_max"]
            )  # maximal yearly generation
            .groupby(n.generators.bus.map(n.buses.country))
            .transform(lambda s: normed(s) * tech_capacities.at[s.name])
            .where(lambda s: s > 0.1, 0.0)
        )  # only capacities above 100kW
        n.generators.loc[tech_i, "p_nom_min"] = n.generators.loc[tech_i, "p_nom"]

        if p_nom_min:
            assert np.isscalar(p_nom_min)
            logger.info(
                f"Scaling capacity stats to {p_nom_min*100:.2f}% of installed capacity acquired from stats."
            )
            n.generators.loc[tech_i, "p_nom_min"] = n.generators.loc[
                tech_i, "p_nom"
            ] * float(p_nom_min)

        if p_nom_max:
            assert np.isscalar(p_nom_max)
            logger.info(
                f"Scaling capacity expansion limit to {p_nom_max*100:.2f}% of installed capacity acquired from stats."
            )
            n.generators.loc[tech_i, "p_nom_max"] = n.generators.loc[
                tech_i, "p_nom_min"
            ] * float(p_nom_max)


def add_nice_carrier_names(n, config):
    carrier_i = n.carriers.index
    nice_names = (
        pd.Series(config["plotting"]["nice_names"])
        .reindex(carrier_i)
        .fillna(carrier_i.to_series().str.title())
    )
    n.carriers["nice_name"] = nice_names
    colors = pd.Series(config["plotting"]["tech_colors"]).reindex(carrier_i)
    if colors.isna().any():
        missing_i = list(colors.index[colors.isna()])
        logger.warning(
            f"tech_colors for carriers {missing_i} not defined " "in config."
        )
    n.carriers["color"] = colors


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("add_electricity")
        sets_path_to_root("pypsa-earth")
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.base_network)
    Nyears = n.snapshot_weightings.objective.sum() / 8760.0

    # Snakemake imports:
    demand_profiles = snakemake.input["demand_profiles"]
    countries = snakemake.config["countries"]
    costs = load_costs(
        tech_costs_path=snakemake.input.tech_costs,
        config_costs=snakemake.config["costs"],
        elec_config=snakemake.config["electricity"],
        countries=countries,
        country_mapping=world_iso,
        Nyears=Nyears,
    ).pipe(rename_technologies, hydrogen_pypsa_mapping)

    ppl = load_powerplants(file_path=snakemake.input.powerplants).pipe(
        map_carriers, carrier_mapping=carrier_pypsa_mapping
    )
    if "renewable_carriers" in snakemake.config["electricity"]:
        renewable_carriers = set(snakemake.config["electricity"]["renewable_carriers"])
    else:
        logger.warning(
            "Missing key `renewable_carriers` under config entry `electricity`. "
            "In future versions, this will raise an error. "
            "Falling back to carriers listed under `renewable`."
        )
        renewable_carriers = set(snakemake.config["renewable"])

    extendable_carriers = snakemake.config["electricity"]["extendable_carriers"]
    if not (set(renewable_carriers) & set(extendable_carriers["Generator"])):
        logger.warning(
            "No renewables found in config entry `extendable_carriers`. "
            "In future versions, these have to be explicitly listed. "
            "Falling back to all renewables."
        )

    conventional_carriers = snakemake.config["electricity"]["conventional_carriers"]
    attach_load(n, demand_profiles)
    update_transmission_costs(
        n, costs, countries, snakemake.config["lines"]["length_factor"]
    )
    conventional_inputs = {
        k: v for k, v in snakemake.input.items() if k.startswith("conventional_")
    }
    attach_conventional_generators(
        n,
        costs,
        ppl,
        conventional_carriers,
        extendable_carriers,
        renewable_carriers,
        snakemake.config.get("conventional", {}),
        conventional_inputs,
        countries,
    )
    attach_wind_and_solar(
        n,
        costs,
        ppl,
        snakemake.input,
        countries,
        renewable_carriers,
        extendable_carriers,
        snakemake.config["lines"]["length_factor"],
    )
    attach_hydro(n, costs, ppl, countries)

    estimate_renewable_capacities_irena(n, snakemake.config)

    update_p_nom_max(n)
    add_nice_carrier_names(n, snakemake.config)

    if not ("weight" in n.generators.columns):
        logger.warning(
            "Unexpected missing 'weight' column, which has been manually added. It may be due to missing generators."
        )
        n.generators["weight"] = pd.Series()

    n.meta = snakemake.config
    n.export_to_netcdf(snakemake.output[0])
