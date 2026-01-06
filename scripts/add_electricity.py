# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Adds electrical generators, load and existing hydro storage units to a base
network.

Relevant Settings
-----------------

.. code:: yaml

    costs:
        year:
        technology_data_version:
        discountrate:
        output_currency:
        country_specific_data:
        cost_scenario:
        financial_case:
        output_currency:
        default_exchange_rate:
        future_exchange_rate_strategy:
        custom_future_exchange_rate:
        rooftop_share:
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
            hydro_max_hours_default:
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

import numpy as np
import pandas as pd
import powerplantmatching as pm
import pypsa
import xarray as xr
from _helpers import (
    apply_currency_conversion,
    build_currency_conversion_cache,
    configure_logging,
    create_logger,
    read_csv_nafix,
    sanitize_carriers,
    sanitize_locations,
    update_p_nom_max,
    get_base_carrier,
)
from powerplantmatching.export import map_country_bus

idx = pd.IndexSlice

logger = create_logger(__name__)


def normed(s):
    return s / s.sum()


def calculate_annuity(n, r):
    """
    Calculate the annuity factor for an asset with lifetime n years and
    discount rate of r, e.g. annuity(20, 0.05) * 20 = 1.6.
    """
    if isinstance(r, pd.Series):
        return pd.Series(1 / n, index=r.index).where(
            r == 0, r / (1.0 - 1.0 / (1.0 + r) ** n)
        )
    elif r > 0:
        return r / (1.0 - 1.0 / (1.0 + r) ** n)
    else:
        return 1 / n


def _add_missing_carriers_from_costs(n, costs, carriers):
    missing_carriers = pd.Index(carriers).difference(n.carriers.index)
    if missing_carriers.empty:
        return

    emissions_cols = (
        costs.columns.to_series().loc[lambda s: s.str.endswith("_emissions")].values
    )
    suptechs = missing_carriers.str.split("-").str[0]
    emissions = costs.loc[suptechs, emissions_cols].fillna(0.0)
    emissions.index = missing_carriers
    n.import_components_from_dataframe(emissions, "Carrier")


def load_costs(tech_costs, config, elec_config, Nyears=1):
    """
    Set all asset costs and other parameters.
    """
    costs = pd.read_csv(tech_costs, index_col=["technology", "parameter"]).sort_index()

    # correct units to MW and output_currency
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3
    costs.unit = costs.unit.str.replace("/kW", "/MW")
    _currency_conversion_cache = build_currency_conversion_cache(
        costs,
        output_currency=config["output_currency"],
        default_exchange_rate=config["default_exchange_rate"],
        future_exchange_rate_strategy=config.get("future_exchange_rate_strategy"),
        custom_future_exchange_rate=config.get("custom_future_rate", None),
    )
    costs = apply_currency_conversion(
        costs,
        config["output_currency"],
        _currency_conversion_cache,
    )

    # apply filter on financial_case and scenario, if they are contained in the cost dataframe
    wished_cost_scenario = config["cost_scenario"]
    wished_financial_case = config["financial_case"]
    for col in ["scenario", "financial_case"]:
        if col in costs.columns:
            costs[col] = costs[col].replace("", pd.NA)

    if "scenario" in costs.columns:
        costs = costs[
            (costs["scenario"].str.casefold() == wished_cost_scenario.casefold())
            | (costs["scenario"].isnull())
        ]

    if "financial_case" in costs.columns:
        costs = costs[
            (costs["financial_case"].str.casefold() == wished_financial_case.casefold())
            | (costs["financial_case"].isnull())
        ]

    costs = costs.value.unstack().fillna(config["fill_values"])

    for attr in ("investment", "lifetime", "FOM", "VOM", "efficiency", "fuel"):
        overwrites = config.get(attr)
        if overwrites is not None:
            overwrites = pd.Series(overwrites)
            costs.loc[overwrites.index, attr] = overwrites
            logger.info(
                f"Overwriting {attr} of {overwrites.index} to {overwrites.values}"
            )

    costs["capital_cost"] = (
        (
            calculate_annuity(costs["lifetime"], costs["discount rate"])
            + costs["FOM"] / 100.0
        )
        * costs["investment"]
        * Nyears
    )

    costs.at["OCGT", "fuel"] = costs.at["gas", "fuel"]
    costs.at["CCGT", "fuel"] = costs.at["gas", "fuel"]

    costs["marginal_cost"] = costs["VOM"] + costs["fuel"] / costs["efficiency"]

    costs = costs.rename(columns={"CO2 intensity": "co2_emissions"})
    # rename because technology data & pypsa earth costs.csv use different names
    # TODO: rename the technologies in hosted tutorial data to match technology data
    costs = costs.rename(
        {
            "hydrogen storage": "hydrogen storage tank",
            "hydrogen storage tank": "hydrogen storage tank",
            "hydrogen storage tank type 1": "hydrogen storage tank",
            "hydrogen underground storage": "hydrogen storage underground",
        },
    )

    costs.at["OCGT", "co2_emissions"] = costs.at["gas", "co2_emissions"]
    costs.at["CCGT", "co2_emissions"] = costs.at["gas", "co2_emissions"]

    costs.at["solar", "capital_cost"] = (
        config["rooftop_share"] * costs.at["solar-rooftop", "capital_cost"]
        + (1 - config["rooftop_share"]) * costs.at["solar-utility", "capital_cost"]
    )
    costs.loc["csp"] = costs.loc["csp-tower"]

    def costs_for_storage(store, link1, link2=None, max_hours=1.0):
        capital_cost = link1["capital_cost"] + max_hours * store["capital_cost"]
        if link2 is not None:
            capital_cost += link2["capital_cost"]
        return pd.Series(
            dict(capital_cost=capital_cost, marginal_cost=0.0, co2_emissions=0.0)
        )

    max_hours = elec_config["max_hours"]
    costs.loc["battery"] = costs_for_storage(
        costs.loc["battery storage"],
        costs.loc["battery inverter"],
        max_hours=max_hours["battery"],
    )
    costs.loc["H2"] = costs_for_storage(
        costs.loc["hydrogen storage tank"],
        costs.loc["fuel cell"],
        costs.loc["electrolysis"],
        max_hours=max_hours["H2"],
    )

    for attr in ("marginal_cost", "capital_cost"):
        overwrites = config.get(attr)
        if overwrites is not None:
            overwrites = pd.Series(overwrites)
            costs.loc[overwrites.index, attr] = overwrites
            logger.info(
                f"Overwriting {attr} of {overwrites.index} to {overwrites.values}"
            )

    return costs


def load_powerplants(
    ppl_fn: str,
    costs: pd.DataFrame = None,
    fill_values: dict = None,
    grouping_years: list = None,
    ) -> pd.DataFrame:
    """
    Load and preprocess powerplant matching data, fill missing datein/dateout, and assign grouping years.
    Parameters
    ----------
    ppl_fn : str
        Path to powerplant matching csv file.
    costs : pd.DataFrame
        DataFrame containing technology costs.
    fill_values : dict
        Dictionary containing default values for lifetime.
    grouping_years : list
        List of years to group build years into.
        
    Returns
    -------
    ppl : pd.DataFrame
        Power plant list DataFrame.
    """
    carrier_dict = {
        "ocgt": "OCGT",
        "ccgt": "CCGT",
        "bioenergy": "biomass",
        "ccgt, thermal": "CCGT",
        "hard coal": "coal",
    }
    ppl = (
        read_csv_nafix(ppl_fn, index_col=0, dtype={"bus": "str"})
        .powerplant.to_pypsa_names()
        .powerplant.convert_country_to_alpha2()
        .rename(columns=str.lower)
        .drop(columns=["efficiency"])
        .replace({"carrier": carrier_dict})
    )
    # drop powerplants with null capacity
    null_ppls = ppl[ppl.p_nom <= 0]
    if not null_ppls.empty:
        logger.warning(f"Drop powerplants with null capacity: {list(null_ppls.name)}.")
        ppl = ppl.drop(null_ppls.index)

    # Fill missing datein and dateout columns
    if costs is not None and fill_values is not None:
        ppl = fill_datein_dateout(ppl, costs, fill_values)

    # Assign grouping years
    if grouping_years is not None:
        ppl["grouping_year"] = get_grouping_year(ppl["datein"], grouping_years)

    return ppl


def fill_datein_dateout(
    ppl: pd.DataFrame, costs: pd.DataFrame, fill_values: dict
) -> pd.DataFrame:
    """
    Fill missing datein and dateout values in ppl DataFrame.

    Parameters
    ----------
    ppl : pd.DataFrame
        Dataframe containing power plants.
    costs : pd.DataFrame
        DataFrame containing cost assumptions.
    fill_values : dict
        Dictionary containing default values for lifetime.

    Returns
    -------
    ppl : pd.DataFrame
        Power plant list DataFrame with filled missing datein and dateout columns.
    """
    # Fill missing datein with mean build year per technology
    if ppl["datein"].isna().any():
        missing_datein = ppl[ppl["datein"].isna()].index
        mean_datein = ppl.groupby("carrier")["datein"].mean()
        ppl.loc[missing_datein, "datein"] = ppl.loc[missing_datein, "carrier"].map(
            mean_datein
        )
        logger.warning(
            f"Filling missing 'datein' for {len(missing_datein)} powerplants with mean build year per technology."
        )

        # Check if there are still missing datein values after carrier mean filling
        if ppl["datein"].isna().any():
            still_missing = ppl[ppl["datein"].isna()].index
            raise ValueError(
                f"Could not fill 'datein' for {len(still_missing)} powerplants by averaging over technology. Please provide explicit 'datein' values for these powerplants."
            )

    # Fill missing dateout based on lifetime from costs DataFrame
    if ppl["dateout"].isna().any():
        missing_dateout = ppl[ppl["dateout"].isna()].index
        ppl.loc[missing_dateout, "dateout"] = ppl.loc[
            missing_dateout, "datein"
        ] + ppl.loc[missing_dateout, "carrier"].map(costs["lifetime"]).fillna(
            fill_values["lifetime"]
        )
        logger.warning(
            f"Filling missing 'dateout' for {len(missing_dateout)} powerplants based on 'datein' and technology lifetimes."
        )

    return ppl


def attach_load(n, demand_profiles):
    """
    Add load profiles to network buses.

    Parameters
    ----------
    n: pypsa network

    demand_profiles: str
        Path to csv file of elecric demand time series, e.g. "resources/demand_profiles.csv"
        Demand profile has snapshots as rows and bus names as columns.

    Returns
    -------
    n : pypsa network
        Now attached with load time series
    """
    demand_df = read_csv_nafix(demand_profiles, index_col=0, parse_dates=True)

    n.madd("Load", demand_df.columns, bus=demand_df.columns, p_set=demand_df)


def attach_dc_costs(lines_or_links, costs, length_factor=1.0, simple_hvdc_costs=False):
    if lines_or_links.empty:
        return

    if lines_or_links.loc[lines_or_links.carrier == "DC"].empty:
        return

    dc_b = lines_or_links.carrier == "DC"
    if simple_hvdc_costs:
        costs = (
            lines_or_links.loc[dc_b, "length"]
            * length_factor
            * costs.at["HVDC overhead", "capital_cost"]
        )
    else:
        costs = (
            lines_or_links.loc[dc_b, "length"]
            * length_factor
            * (
                (1.0 - lines_or_links.loc[dc_b, "underwater_fraction"])
                * costs.at["HVDC overhead", "capital_cost"]
                + lines_or_links.loc[dc_b, "underwater_fraction"]
                * costs.at["HVDC submarine", "capital_cost"]
            )
            + costs.at["HVDC inverter pair", "capital_cost"]
        )
    lines_or_links.loc[dc_b, "capital_cost"] = costs


def update_transmission_costs(n, costs, length_factor=1.0, simple_hvdc_costs=False):
    n.lines["capital_cost"] = (
        n.lines["length"] * length_factor * costs.at["HVAC overhead", "capital_cost"]
    )

    attach_dc_costs(
        lines_or_links=n.links,
        costs=costs,
        length_factor=length_factor,
        simple_hvdc_costs=simple_hvdc_costs,
    )
    attach_dc_costs(
        lines_or_links=n.lines,
        costs=costs,
        length_factor=length_factor,
        simple_hvdc_costs=simple_hvdc_costs,
    )


def get_grouping_year(build_year, grouping_years):
    """
    Map build_year to the nearest grouping year (rounded up).
    
    Example:
        grouping_years = [1980, 2000, 2010, 2015, 2020]
        build_year = 2012 → returns 2015
        build_year = 2018 → returns 2020
    """
    indices = np.digitize(build_year, grouping_years, right=True)
    return np.take(grouping_years, indices)


def aggregate_ppl_by_bus_carrier_year(ppl: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate power plants by (bus, carrier, grouping_year).
    
    Creates a new carrier name with grouping year suffix (e.g., "CCGT-2020")
    and aggregates capacity and other attributes.
    
    Parameters
    ----------
    ppl : pd.DataFrame
        Power plant DataFrame with columns: bus, carrier, grouping_year, 
        p_nom, efficiency, marginal_cost, datein, dateout and so on.
    
    Returns
    -------
    pd.DataFrame
        Aggregated power plants with columns: bus, carrier, carrier_gy,
        p_nom, efficiency, marginal_cost, build_year, lifetime.
        
    Example
    -------
    Input:
        bus    carrier  grouping_year  p_nom
        bus1   CCGT     2015           100
        bus1   CCGT     2015           200
        bus1   CCGT     2020           150
        
    Output:
        bus    carrier  carrier_gy   p_nom
        bus1   CCGT     CCGT-2015    300
        bus1   CCGT     CCGT-2020    150
    """
    # Add grouping year to carrier name
    ppl = ppl.copy()
    ppl["carrier_gy"] = ppl["carrier"] + "-" + ppl["grouping_year"].astype(str)
    
    # Group by (bus, carrier_gy) and aggregate
    agg_dict = {
        "carrier": "first",
        "p_nom": "sum",
        "datein": "mean",
        "dateout": "mean",
        "max_hours": "mean",
    }

    if "marginal_cost" in ppl.columns:
        agg_dict["marginal_cost"] = "mean"
    if "capital_cost" in ppl.columns:
        agg_dict["capital_cost"] = "mean"
    if "efficiency" in ppl.columns:
        agg_dict["efficiency"] = "mean"
    if "country" in ppl.columns:
        agg_dict["country"] = "first"
    
    ppl_grouped = ppl.groupby(["bus", "carrier_gy"]).agg(agg_dict).reset_index()
    
    # Calculate build_year and lifetime
    ppl_grouped["build_year"] = ppl_grouped["datein"].astype(int)
    ppl_grouped["lifetime"] = (ppl_grouped["dateout"] - ppl_grouped["datein"]).fillna(np.inf)
    
    # Set index as "bus carrier_gy"
    ppl_grouped = ppl_grouped.set_index(ppl_grouped["bus"] + " " + ppl_grouped["carrier_gy"])

    return ppl_grouped


def attach_wind_and_solar(
    n: pypsa.Network,
    costs: pd.DataFrame,
    ppl: pd.DataFrame,
    input_files: dict,  # snakemake input
    technologies: set,
    extendable_carriers: dict,
    line_length_factor: float,
):
    """
    Add existing and extendable wind and solar generators to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to modify.
    costs : pd.DataFrame
        DataFrame containing technology costs.
    ppl : pd.DataFrame
        Power plant DataFrame.
    input_files : dict
        Snakemake input object.
    technologies : set
        Set of renewable technologies to be added.
    extendable_carriers : dict
        Dictionary of extendable carriers for different component types.
    line_length_factor : float
        Factor to adjust line lengths for connection cost calculations.

    Returns
    -------
    None
    """
    # TODO: rename tech -> carrier, technologies -> carriers
    _add_missing_carriers_from_costs(n, costs, technologies)

    df = ppl.rename(columns={"country": "Country"})

    for tech in technologies:
        if tech == "hydro":
            continue

        if tech == "offwind-ac":
            # add all offwind wind power plants by default as offwind-ac
            df["carrier"] = df["carrier"].mask(
                df.technology == "Offshore", "offwind-ac"
            )

        df["carrier"] = df["carrier"].mask(df.technology == "Onshore", "onwind")

        with xr.open_dataset(getattr(input_files, "profile_" + tech)) as ds:
            if ds.indexes["bus"].empty:
                continue

            suptech = tech.split("-", 2)[0]
            if suptech == "offwind":
                underwater_fraction = ds["underwater_fraction"].to_pandas()
                connection_cost = (
                    line_length_factor
                    * ds["average_distance"].to_pandas()
                    * (
                        underwater_fraction
                        * costs.at[tech + "-connection-submarine", "capital_cost"]
                        + (1.0 - underwater_fraction)
                        * costs.at[tech + "-connection-underground", "capital_cost"]
                    )
                )
                capital_cost = (
                    costs.at["offwind", "capital_cost"]
                    + costs.at[tech + "-station", "capital_cost"]
                    + connection_cost
                )
                logger.info(
                    "Added connection cost of {:0.0f}-{:0.0f} {}/MW/a to {}".format(
                        connection_cost.min(),
                        connection_cost.max(),
                        snakemake.params.costs["output_currency"],
                        tech,
                    )
                )
            else:
                capital_cost = costs.at[tech, "capital_cost"]

            # Get p_max_pu profile for all buses
            p_max_pu = ds["profile"].transpose("time", "bus").to_pandas()
            p_nom_max_per_bus = ds["p_nom_max"].to_pandas()
            weight = ds["weight"].to_pandas()

            # Add existing generators (not extendable)
            tech_ppl = df.query("carrier == @tech")
            if not tech_ppl.empty:
                buses = n.buses.loc[ds.indexes["bus"]]
                existing = map_country_bus(tech_ppl, buses)

                # Aggregate existing generators by (bus, carrier, grouping_year)
                existing_grouped = aggregate_ppl_by_bus_carrier_year(existing)

                # Add carrier with grouping year
                for carrier_gy in existing_grouped["carrier_gy"].unique():
                    if carrier_gy not in n.carriers.index:
                        n.add(
                            "Carrier",
                            carrier_gy,
                            co2_emissions=n.carriers.at[tech, "co2_emissions"]
                            if tech in n.carriers.index else 0.0,
                        )

                # Get p_max_pu for each existing generator's bus
                existing_p_max_pu = p_max_pu[existing_grouped["bus"].values]
                existing_p_max_pu.columns = existing_grouped.index

                n.madd(
                    "Generator",
                    existing_grouped.index,
                    bus=existing_grouped["bus"],
                    carrier=existing_grouped["carrier_gy"],
                    p_nom=existing_grouped["p_nom"],
                    p_nom_extendable=False,
                    p_nom_min=existing_grouped["p_nom"],
                    p_nom_max=existing_grouped["p_nom"],
                    p_max_pu=existing_p_max_pu,
                    marginal_cost=costs.at[suptech, "marginal_cost"],
                    capital_cost=capital_cost,
                    efficiency=costs.at[suptech, "efficiency"],
                    build_year=existing_grouped["build_year"],
                    lifetime=existing_grouped["lifetime"],
                )

                # Calculate existing capacity per bus for adjusting p_nom_max
                existing_per_bus = existing.groupby("bus")["p_nom"].sum()
                existing_per_bus = existing_per_bus.reindex(ds.indexes["bus"]).fillna(0)

                logger.info(
                    f"Added {len(existing)} existing {tech} generators "
                    f"with total capacity {existing['p_nom'].sum() / 1e3:.2f} GW"
                )
            else:
                existing_per_bus = pd.Series(0.0, index=ds.indexes["bus"])

            # Add extendable generators
            if tech in extendable_carriers["Generator"]:
                # Adjust p_nom_max by subtracting existing capacities
                adjusted_p_nom_max = (p_nom_max_per_bus - existing_per_bus).clip(
                    lower=0.0
                )

                n.madd(
                    "Generator",
                    ds.indexes["bus"],
                    suffix=" " + tech,
                    bus=ds.indexes["bus"],
                    carrier=tech,
                    p_nom_extendable=True,
                    p_nom_max=adjusted_p_nom_max,
                    p_max_pu=p_max_pu,
                    weight=weight,
                    marginal_cost=costs.at[suptech, "marginal_cost"],
                    capital_cost=capital_cost,
                    efficiency=costs.at[suptech, "efficiency"],
                    lifetime=costs.at[suptech, "lifetime"],
                )

                logger.info(
                    f"Added extendable {tech} generators at {len(ds.indexes['bus'])} buses "
                    f"with total potential {adjusted_p_nom_max.sum() / 1e3:.2f} GW"
                )


def attach_conventional_generators(
    n: pypsa.Network,
    costs: pd.DataFrame,
    ppl: pd.DataFrame,
    conventional_carriers: list,
    extendable_carriers: dict,
    renewable_carriers: set,
    conventional_config: list,
    conventional_inputs: dict,
):
    """
    Add existing conventional generators to the network and extendable conventional generators at all buses.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to modify.
    costs : pd.DataFrame
        DataFrame containing technology costs.
    ppl : pd.DataFrame
        Power plant DataFrame.
    conventional_carriers : list
        List of conventional carriers to be added.
    extendable_carriers : dict
        Dictionary of extendable carriers for different component types.
    renewable_carriers : set
        Set of renewable carriers.
    conventional_config : list
        List of conventional configuration settings.
    conventional_inputs : dict
        Dictionary of conventional input parameters.

    Returns
    -------
    None
    """
    carriers = set(conventional_carriers) | (
        set(extendable_carriers["Generator"]) - set(renewable_carriers)
    )
    _add_missing_carriers_from_costs(n, costs, carriers)

    ppl = (
        ppl.query("carrier in @carriers")
        .join(costs, on="carrier", rsuffix="_r")
        .rename(index=lambda s: "C" + str(s))
    )
    ppl["efficiency"] = ppl.efficiency.fillna(ppl.efficiency)

    # Aggregate power plants by (bus, carrier, grouping_year)
    ppl_grouped = aggregate_ppl_by_bus_carrier_year(ppl)

    # Add carriers with grouping year
    for carrier_gy in ppl_grouped["carrier_gy"].unique():
        base_carrier = carrier_gy.rsplit("-", 1)[0]
        if carrier_gy not in n.carriers.index:
            n.add(
                "Carrier",
                carrier_gy,
                co2_emissions=n.carriers.at[base_carrier, "co2_emissions"],
            )
    
    logger.info(
        "Adding {} existing generators with capacities [GW] \n{}".format(
            len(ppl), ppl.groupby("carrier").p_nom.sum().div(1e3).round(2)
        )
    )

    n.madd(
        "Generator",
        ppl_grouped["bus"] + " " + ppl_grouped["carrier_gy"],
        carrier=ppl_grouped["carrier_gy"],
        bus=ppl_grouped["bus"],
        p_nom=ppl_grouped["p_nom"],
        p_nom_extendable=False,
        p_nom_min=ppl_grouped["p_nom"],
        p_nom_max=ppl_grouped["p_nom"],
        efficiency=ppl_grouped["efficiency"],
        marginal_cost=ppl_grouped["marginal_cost"],
        capital_cost=ppl_grouped["capital_cost"],
        build_year=ppl_grouped["build_year"],
        lifetime=ppl_grouped["lifetime"],
    )

    # Add extendable conventional generators
    extendable_conventional = set(extendable_carriers["Generator"]) - set(renewable_carriers)

    for carrier in extendable_conventional:
        carrier_buses = ppl[ppl.carrier == carrier]["bus"].unique()
        n.madd(
            "Generator",
            carrier_buses,
            suffix=" " + carrier,
            carrier=carrier,
            bus=carrier_buses,
            p_nom_extendable=True,
            efficiency=costs.at[carrier, "efficiency"],
            marginal_cost=costs.at[carrier, "marginal_cost"],
            capital_cost=costs.at[carrier, "capital_cost"],
            lifetime=costs.at[carrier, "lifetime"],
        )

    logger.info(
        f"Added extendable {extendable_conventional} generators at {len(carrier_buses)} buses."
    )

    for carrier in conventional_config:
        # Generators with technology affected
        idx = n.generators.query("carrier == @carrier").index

        for attr in list(set(conventional_config[carrier]) & set(n.generators)):
            values = conventional_config[carrier][attr]

            if f"conventional_{carrier}_{attr}" in conventional_inputs:
                # Values affecting generators of technology k country-specific
                # First map generator buses to countries; then map countries to p_max_pu
                values = read_csv_nafix(values, index_col=0).iloc[:, 0]
                bus_values = n.buses.country.map(values)
                n.generators[attr].update(
                    n.generators.loc[idx].bus.map(bus_values).dropna()
                )
            else:
                # Single value affecting all generators of technology k indiscriminantely of country
                n.generators.loc[idx, attr] = values


def attach_hydro(n: pypsa.Network, costs: pd.DataFrame, ppl: pd.DataFrame) -> None:
    """
    Add existing hydro powerplants to the network as Hydro Storage units, Run-Of-River generators, and Pumped Hydro storage units.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to modify.
    costs : pd.DataFrame
        DataFrame containing technology costs.
    ppl : pd.DataFrame
        Power plant DataFrame.

    Returns
    -------
    None
    """
    if "hydro" not in snakemake.params.renewable:
        return
    c = snakemake.params.renewable["hydro"]
    carriers = c.get("carriers", ["ror", "PHS", "hydro"])

    _add_missing_carriers_from_costs(n, costs, carriers)

    ppl = (
        ppl.query('carrier == "hydro"')
        .assign(ppl_id=lambda df: df.index)
        .reset_index(drop=True)
        .rename(index=lambda s: str(s) + " hydro")
    )

    # Current fix, NaN technologies set to ROR
    if ppl.technology.isna().any():
        n_nans = ppl.technology.isna().sum()
        logger.warning(
            f"Identified {n_nans} hydro powerplants with unknown technology.\n"
            "Initialized to 'Run-Of-River'"
        )
        ppl.loc[ppl.technology.isna(), "technology"] = "Run-Of-River"

    # Map technology to carrier before aggregation
    tech_to_carrier = {
        "Run-Of-River": "ror",
        "Pumped Storage": "PHS",
        "Reservoir": "hydro",
    }
    ppl["carrier"] = ppl["technology"].map(tech_to_carrier)

    # Aggregate by (bus, carrier, grouping_year)
    ppl_grouped = aggregate_ppl_by_bus_carrier_year(ppl)

    ror = ppl_grouped[ppl_grouped["carrier"] == "ror"]
    phs = ppl_grouped[ppl_grouped["carrier"] == "PHS"]
    hydro = ppl_grouped[ppl_grouped["carrier"] == "hydro"]

    inflow_idx = ror.index.union(hydro.index)
    if not inflow_idx.empty:
        with xr.open_dataarray(snakemake.input.profile_hydro) as inflow:
            found_plants = ppl.ppl_id[ppl.ppl_id.isin(inflow.indexes["plant"])]
            missing_plants_idxs = ppl.index.difference(found_plants.index)

            # if missing time series are found, notify the user and exclude missing hydro plants
            if not missing_plants_idxs.empty:
                # original total p_nom
                total_p_nom = ror.p_nom.sum() + hydro.p_nom.sum()

                ror = ror.loc[ror.index.intersection(found_plants.index)]
                hydro = hydro.loc[hydro.index.intersection(found_plants.index)]
                # loss of p_nom
                loss_p_nom = ror.p_nom.sum() + hydro.p_nom.sum() - total_p_nom

                logger.warning(
                    f"'{snakemake.input.profile_hydro}' is missing inflow time-series for at least one bus: {', '.join(missing_plants_idxs)}."
                    f"Corresponding hydro plants are dropped, corresponding to a total loss of {loss_p_nom:.2f}MW out of {total_p_nom:.2f}MW."
                )

            # if there are any plants for which runoff data are available
            if not found_plants.empty:
                inflow_t = (
                    inflow.sel(plant=found_plants.values)
                    .assign_coords(plant=found_plants.index)
                    .rename({"plant": "name"})
                    .transpose("time", "name")
                    .to_pandas()
                )

                # Aggregate inflow by (bus, carrier, grouping_year)
                inflow_dict = {}
                for idx in ppl_grouped.index:
                    bus = ppl_grouped.at[idx, "bus"]
                    carrier_gy = ppl_grouped.at[idx, "carrier_gy"]
                    grouping_year = int(carrier_gy.rsplit("-", 1)[1])
                    carrier = ppl_grouped.at[idx, "carrier"]

                    # Find original plants in this group
                    mask = (
                        (ppl["bus"] == bus) &
                        (ppl["carrier"] == carrier) &
                        (ppl["grouping_year"] == grouping_year)
                    )
                    original_plants = ppl[mask].index
                    valid_plants = original_plants[original_plants.isin(inflow_t.columns)]

                    if not valid_plants.empty:
                        inflow_dict[idx] = inflow_t[valid_plants].sum(axis=1)

                inflow_agg = pd.DataFrame(inflow_dict, index=inflow_t.index)

    # Add carriers with grouping year
    for carrier_gy in ppl_grouped["carrier_gy"].unique():
        if carrier_gy not in n.carriers.index:
            n.add("Carrier", carrier_gy, co2_emissions=0.0)

    if "ror" in carriers and not ror.empty:
        n.madd(
            "Generator",
            ror.index,
            carrier=ror["carrier_gy"],
            bus=ror["bus"],
            p_nom=ror["p_nom"],
            p_nom_extendable=False,
            efficiency=costs.at["ror", "efficiency"],
            capital_cost=costs.at["ror", "capital_cost"],
            weight=ror["p_nom"],
            p_max_pu=(
                inflow_agg[ror.index]
                .divide(ror["p_nom"], axis=1)
                .where(lambda df: df <= 1.0, other=1.0)
            ),
            build_year=ror["build_year"],
            lifetime=ror["lifetime"],
        )

        logger.info(
            f"Added {len(ror)} ror generators with {ror['p_nom'].sum() / 1e3:.2f} GW"
        )

    if "PHS" in carriers and not phs.empty:
        # fill missing max hours to config value and
        # assume no natural inflow due to lack of data
        phs = phs.replace({"max_hours": {0: c["PHS_max_hours"]}})
        n.madd(
            "StorageUnit",
            phs.index,
            carrier=phs["carrier_gy"],
            bus=phs["bus"],
            p_nom=phs["p_nom"],
            p_nom_extendable=False,
            capital_cost=costs.at["PHS", "capital_cost"],
            max_hours=phs["max_hours"],
            efficiency_store=np.sqrt(costs.at["PHS", "efficiency"]),
            efficiency_dispatch=np.sqrt(costs.at["PHS", "efficiency"]),
            cyclic_state_of_charge=True,
            build_year=phs["build_year"],
            lifetime=phs["lifetime"],
        )

        logger.info(
            f"Added {len(phs)} PHS storage units with {phs['p_nom'].sum() / 1e3:.2f} GW"
        )

    if "hydro" in carriers and not hydro.empty:
        hydro_max_hours = c.get("hydro_max_hours")
        hydro_stats = (
            pd.read_csv(
                snakemake.input.hydro_capacities,
                comment="#",
                na_values=["-"],
                index_col=0,
            )
            .groupby("Country")
            .sum()
        )
        e_target = hydro_stats["E_store[TWh]"].clip(lower=0.2) * 1e6
        e_installed = hydro.eval("p_nom * max_hours").groupby(hydro.country).sum()
        e_missing = e_target - e_installed
        missing_mh_i = hydro.query("max_hours.isnull()").index

        if hydro_max_hours == "energy_capacity_totals_by_country":
            max_hours_country = (
                e_missing / hydro.loc[missing_mh_i].groupby("country").p_nom.sum()
            )

        elif hydro_max_hours == "estimate_by_large_installations":
            max_hours_country = (
                hydro_stats["E_store[TWh]"] * 1e3 / hydro_stats["p_nom_discharge[GW]"]
            )

        max_hours_country.clip(lower=0, inplace=True)

        missing_countries = pd.Index(hydro["country"].unique()).difference(
            max_hours_country.dropna().index
        )
        if not missing_countries.empty:
            logger.warning(
                "Assuming max_hours=6 for hydro reservoirs in the countries: {}".format(
                    ", ".join(missing_countries)
                )
            )
        hydro_max_hours_default = c.get("hydro_max_hours_default", 6.0)
        hydro_max_hours = hydro.max_hours.where(
            hydro.max_hours > 0, hydro.country.map(max_hours_country)
        ).fillna(hydro_max_hours_default)

        n.madd(
            "StorageUnit",
            hydro.index,
            carrier=hydro["carrier_gy"],
            bus=hydro["bus"],
            p_nom=hydro["p_nom"],
            p_nom_extendable=False,
            max_hours=hydro_max_hours,
            capital_cost=(
                costs.at["hydro", "capital_cost"]
                if c.get("hydro_capital_cost")
                else 0.0
            ),
            marginal_cost=costs.at["hydro", "marginal_cost"],
            p_max_pu=1.0,  # dispatch
            p_min_pu=0.0,  # store
            efficiency_dispatch=costs.at["hydro", "efficiency"],
            efficiency_store=0.0,
            cyclic_state_of_charge=True,
            inflow=inflow_agg[hydro.index],
            build_year=hydro["build_year"],
            lifetime=hydro["lifetime"],
        )

        logger.info(
            f"Added {len(hydro)} hydro storage units with {hydro['p_nom'].sum() / 1e3:.2f} GW"
        )


def attach_extendable_generators(
    n: pypsa.Network, costs: pd.DataFrame, ppl: pd.DataFrame
) -> None:
    """
    Add extendable conventional generators (OCGT, CCGT, nuclear) with zero capacity.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to modify.
    costs : pd.DataFrame
        DataFrame containing technology costs.
    ppl : pd.DataFrame
        Power plant DataFrame.

    Returns
    -------
    None
    """
    logger.warning("The function is deprecated with the next release")
    elec_opts = snakemake.params.electricity
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


def estimate_renewable_capacities_irena(
    n: pypsa.Network, estimate_renewable_capacities_config: dict, countries_config: list
) -> None:
    """
    Estimate renewable capacities based on IRENA statistics.

    This function distributes country-level IRENA capacity statistics to
    extendable generators, accounting for existing capacity already in the network.
    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to modify.
    estimate_renewable_capacities_config : dict
        Configuration dictionary for renewable capacity estimation.
    countries_config : list
        List of countries to consider for capacity estimation.

    Returns
    -------
    None
    """
    stats = estimate_renewable_capacities_config["stats"]
    if not stats:
        return

    year = estimate_renewable_capacities_config["year"]
    tech_map = estimate_renewable_capacities_config["technology_mapping"]
    tech_keys = list(tech_map.keys())
    countries = countries_config

    p_nom_max = estimate_renewable_capacities_config["p_nom_max"]
    p_nom_min = estimate_renewable_capacities_config["p_nom_min"]

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

        # Get base carrier for each generator
        base_carriers = n.generators["carrier"].apply(get_base_carrier)

        # Separate existing and extendable generators
        existing_i = n.generators[
            base_carriers.isin(techs) & (n.generators["carrier"] != base_carriers)
        ].index
        extendable_i = n.generators[n.generators["carrier"].isin(techs)].index

        # Calculate existing capacities per country
        existing_capacity_per_country = (
            n.generators.loc[existing_i, "p_nom"]
            .groupby(n.generators.loc[existing_i, "bus"].map(n.buses.country))
            .sum()
            .reindex(countries, fill_value=0.0)
        )

        # Calculate remaining capacity to distribute (IRENA - existing)
        remaining_capacities = (tech_capacities - existing_capacity_per_country).clip(
            lower=0.0
        )

        logger.info(
            f"For {ppm_technology}: {stats} total = {tech_capacities.sum():.1f} MW, "
            f"Existing = {existing_capacity_per_country.sum():.1f} MW, "
            f"Remaining to distribute = {remaining_capacities.sum():.1f} MW"
        )

        # Skip if no remaining capacity to distribute
        if remaining_capacities.sum() <= 0.0:
            logger.info(
                f"Existing capacity exceeds or equals {stats} stats for {ppm_technology}. "
                "No additional capacity distributed to extendable generators."
            )
            continue

        # Skip if no extendable generators for the technology
        if extendable_i.empty:
            logger.info(
                f"No extendable generators for {ppm_technology}. "
                "Cannot distribute remaining capacity."
            )
            continue

        # Distribute remaining capacity to extendable generators based on potential
        potential = (
            n.generators_t.p_max_pu[extendable_i].mean()
            * n.generators.loc[extendable_i, "p_nom_max"]
        )

        # Distribute by country, weighted by potential
        distributed_capacities = (
            potential.groupby(
                n.generators.loc[extendable_i, "bus"].map(n.buses.country)
            )
            .transform(
                lambda s: (
                    normed(s) * remaining_capacities.at[s.name]
                    if s.name in remaining_capacities.index
                    else 0.0
                )
            )
            .where(lambda s: s > 0.1, 0.0)
        )  # only capacities above 100kW

        n.generators.loc[extendable_i, "p_nom"] = distributed_capacities
        n.generators.loc[extendable_i, "p_nom_min"] = distributed_capacities

        if p_nom_min:
            assert np.isscalar(p_nom_min)
            logger.info(
                f"Scaling capacity stats to {p_nom_min*100:.2f}% of installed capacity acquired from stats."
            )
            n.generators.loc[extendable_i, "p_nom_min"] = n.generators.loc[
                extendable_i, "p_nom"
            ] * float(p_nom_min)

        if p_nom_max:
            assert np.isscalar(p_nom_max)
            logger.info(
                f"Scaling capacity expansion limit to {p_nom_max*100:.2f}% of installed capacity acquired from stats."
            )
            n.generators.loc[extendable_i, "p_nom_max"] = n.generators.loc[
                extendable_i, "p_nom_min"
            ] * float(p_nom_max)


def add_nice_carrier_names(n, config):
    """
    Add nice names and colors to carriers.

    For vintage carriers (e.g., "solar-2020"), uses the nice name and color
    from the base carrier (e.g., "solar").
    """
    carrier_i = n.carriers.index

    # Get base carriers (handles "solar-2020" -> "solar")
    base_carriers = carrier_i.to_series().apply(get_base_carrier)

    # Map nice names from base carrier
    nice_names_config = pd.Series(config["plotting"]["nice_names"])
    nice_names = (
        base_carriers.map(nice_names_config)
        .fillna(carrier_i.to_series().str.title())
    )
    n.carriers["nice_name"] = nice_names

    # Map colors from base carrier
    colors_config = pd.Series(config["plotting"]["tech_colors"])
    colors = base_carriers.map(colors_config)

    if colors.isna().any():
        missing_i = list(colors.index[colors.isna()])
        logger.warning(
            f"tech_colors for carriers {missing_i} not defined " "in config."
        )
    n.carriers["color"] = colors


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("add_electricity")

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.base_network)
    Nyears = n.snapshot_weightings.objective.sum() / 8760.0

    # Snakemake imports:
    demand_profiles = snakemake.input["demand_profiles"]

    costs = load_costs(
        snakemake.input.tech_costs,
        snakemake.params.costs,
        snakemake.params.electricity,
        Nyears,
    )
    ppl = load_powerplants(
        snakemake.input.powerplants,
        costs,
        snakemake.params.costs["fill_values"],
        snakemake.params.existing_capacities["grouping_years_power"],
    )

    if "renewable_carriers" in snakemake.params.electricity:
        renewable_carriers = set(snakemake.params.electricity["renewable_carriers"])
    else:
        logger.warning(
            "Missing key `renewable_carriers` under config entry `electricity`. "
            "In future versions, this will raise an error. "
            "Falling back to carriers listed under `renewable`."
        )
        renewable_carriers = set(snakemake.params.renewable)

    extendable_carriers = snakemake.params.electricity["extendable_carriers"]
    if not (set(renewable_carriers) & set(extendable_carriers["Generator"])):
        logger.warning(
            "No renewables found in config entry `extendable_carriers`. "
            "In future versions, these have to be explicitly listed. "
            "Falling back to all renewables."
        )

    conventional_carriers = snakemake.params.electricity["conventional_carriers"]
    attach_load(n, demand_profiles)
    update_transmission_costs(n, costs, snakemake.params.length_factor)
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
        snakemake.params.conventional,
        conventional_inputs,
    )
    attach_wind_and_solar(
        n,
        costs,
        ppl,
        snakemake.input,
        renewable_carriers,
        extendable_carriers,
        snakemake.params.length_factor,
    )
    attach_hydro(n, costs, ppl)

    if snakemake.params.electricity.get("estimate_renewable_capacities"):
        estimate_renewable_capacities_irena(
            n,
            snakemake.params.electricity["estimate_renewable_capacities"],
            snakemake.params.countries,
        )

    update_p_nom_max(n)
    add_nice_carrier_names(n, snakemake.config)

    if not ("weight" in n.generators.columns):
        logger.warning(
            "Unexpected missing 'weight' column, which has been manually added. It may be due to missing generators."
        )
        n.generators["weight"] = pd.Series()

    sanitize_carriers(n, snakemake.config)
    if "location" in n.buses:
        sanitize_locations(n)

    n.meta = snakemake.config
    n.export_to_netcdf(snakemake.output[0])
