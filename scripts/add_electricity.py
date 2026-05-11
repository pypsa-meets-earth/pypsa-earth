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
        output_currency:

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
    configure_logging,
    create_logger,
    read_csv_nafix,
    sanitize_carriers,
    sanitize_locations,
    update_p_nom_max,
)

idx = pd.IndexSlice

logger = create_logger(__name__)


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


def get_grouping_year(build_year: int, grouping_years: list) -> int:
    """
    Map build_year to the nearest grouping year bin.

    Each build year is assigned to the first grouping year that is greater
    than or equal to it (i.e. the right edge of its bin).

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
    ppl_grouped["lifetime"] = (ppl_grouped["dateout"] - ppl_grouped["datein"]).fillna(
        np.inf
    )

    # Set index as "bus carrier_gy"
    ppl_grouped = ppl_grouped.set_index(
        ppl_grouped["bus"] + " " + ppl_grouped["carrier_gy"]
    )

    return ppl_grouped


def aggregate_inflow_by_group(
    ppl: pd.DataFrame,
    ppl_grouped: pd.DataFrame,
    inflow_t: pd.DataFrame,
) -> pd.DataFrame:
    """
    Aggregate inflow time series by (bus, carrier, grouping_year) groups.

    Parameters
    ----------
    ppl : pd.DataFrame
        Original (ungrouped) power plant DataFrame with columns: bus, carrier, grouping_year.
    ppl_grouped : pd.DataFrame
        Aggregated power plant DataFrame with columns: bus, carrier, carrier_gy.
    inflow_t : pd.DataFrame
        Inflow time series DataFrame with plant indices as columns.

    Returns
    -------
    pd.DataFrame
        Aggregated inflow time series with ppl_grouped indices as columns.
    """
    inflow_dict = {}
    for idx in ppl_grouped.index:
        bus = ppl_grouped.at[idx, "bus"]
        carrier_gy = ppl_grouped.at[idx, "carrier_gy"]
        grouping_year = int(carrier_gy.rsplit("-", 1)[1])
        carrier = ppl_grouped.at[idx, "carrier"]

        mask = (
            (ppl["bus"] == bus)
            & (ppl["carrier"] == carrier)
            & (ppl["grouping_year"] == grouping_year)
        )
        original_plants = ppl[mask].index
        valid_plants = original_plants[original_plants.isin(inflow_t.columns)]

        if not valid_plants.empty:
            inflow_dict[idx] = inflow_t[valid_plants].sum(axis=1)

    return pd.DataFrame(inflow_dict, index=inflow_t.index)


def get_irena_targets_for_carrier(
    carrier: str,
    estimate_renewable_capacities_config: dict,
    countries: list,
) -> pd.Series:
    """
    Return IRENA installed capacity targets for a given carrier as a Series
    indexed by country (MW).

    The function reads IRENASTAT installed capacity data using the existing
    `estimate_renewable_capacities` configuration.
    Offshore wind is mapped entirely to offwind-ac.
    """
    if not estimate_renewable_capacities_config:
        return pd.Series(dtype=float)

    stats = estimate_renewable_capacities_config.get("stats")
    year = estimate_renewable_capacities_config.get("year")
    tech_map = estimate_renewable_capacities_config.get("technology_mapping", {})

    if not stats or stats != "irena":
        return pd.Series(dtype=float)

    if not countries or not tech_map:
        return pd.Series(dtype=float)

    # Read IRENASTAT exactly as in the legacy implementation
    capacities = pm.data.IRENASTAT().powerplant.convert_country_to_alpha2()

    missing = list(set(countries).difference(capacities.Country.unique()))
    if missing:
        logger.info(
            f"The countries {missing} are not provided in the IRENA stats and hence not scaled"
        )

    capacities = capacities.query(
        "Year == @year and Technology in @tech_map.keys() and Country in @countries"
    )

    if capacities.empty:
        return pd.Series(0.0, index=countries)

    capacities = capacities.groupby(["Technology", "Country"]).Capacity.sum()

    targets = pd.Series(0.0, index=countries, dtype=float)

    for ppm_technology, techs in tech_map.items():
        # Offshore is entirely mapped to offwind-ac
        mapped_techs = [
            "offwind-ac" if tech == "offwind-dc" else tech for tech in techs
        ]

        if carrier not in mapped_techs:
            continue

        if ppm_technology not in capacities.index.get_level_values(0):
            logger.info(
                f"technology {ppm_technology} is not provided by {stats} and therefore not estimated"
            )
            continue

        targets = targets.add(
            capacities.loc[ppm_technology].reindex(countries, fill_value=0.0),
            fill_value=0.0,
        )

    return targets


def attach_wind_and_solar(
    n: pypsa.Network,
    costs: pd.DataFrame,
    ppl: pd.DataFrame,
    input_files: dict,
    carriers: set,
    extendable_carriers: dict,
    line_length_factor: float,
) -> None:
    """
    Attach wind and solar generators.

    Existing capacities are taken from powerplants.csv and spatialized to buses.
    National capacity gaps with respect to IRENA targets are filled and
    redistributed uniformly across buses within each country.

    Offshore wind is treated entirely as offwind-ac.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to modify.
    costs : pd.DataFrame
        DataFrame containing technology costs.
    ppl : pd.DataFrame
        Power plant DataFrame.
    input_files : dict
        Snakemake input object containing renewable profile files.
    carriers : set
        Set of renewable carriers to be added.
    extendable_carriers : dict
        Dictionary of extendable carriers for different component types.
    line_length_factor : float
        Factor to adjust line lengths for connection cost calculations.

    Returns
    -------
    None
    """
    _add_missing_carriers_from_costs(n, costs, carriers)

    df = ppl.rename(columns={"country": "Country"}).copy()

    # Force offshore to offwind-ac
    df.loc[df.technology == "Offshore", "carrier"] = "offwind-ac"
    df.loc[df.technology == "Onshore", "carrier"] = "onwind"

    estimate_cfg = snakemake.params.electricity.get("estimate_renewable_capacities", {})
    countries = snakemake.params.countries

    for carrier in carriers:
        if carrier == "hydro":
            continue

        with xr.open_dataset(getattr(input_files, "profile_" + carrier)) as ds:
            if ds.indexes["bus"].empty:
                continue

            supcarrier = carrier.split("-", 2)[0]

            if supcarrier == "offwind":
                underwater_fraction = ds["underwater_fraction"].to_pandas()
                connection_cost = (
                    line_length_factor
                    * ds["average_distance"].to_pandas()
                    * (
                        underwater_fraction
                        * costs.at[carrier + "-connection-submarine", "capital_cost"]
                        + (1.0 - underwater_fraction)
                        * costs.at[carrier + "-connection-underground", "capital_cost"]
                    )
                )
                capital_cost = (
                    costs.at["offwind", "capital_cost"]
                    + costs.at[carrier + "-station", "capital_cost"]
                    + connection_cost
                )

                logger.info(
                    "Added connection cost of {:0.0f}-{:0.0f} {}/MW/a to {}".format(
                        connection_cost.min(),
                        connection_cost.max(),
                        snakemake.config["costs"]["output_currency"],
                        carrier,
                    )
                )
            else:
                capital_cost = costs.at[carrier, "capital_cost"]

            # Existing capacity from powerplants.csv
            ppl_carrier = df.query("carrier == @carrier")

            if not ppl_carrier.empty:
                valid = ppl_carrier[ppl_carrier.bus.isin(n.buses.index)].copy()
                missing = ppl_carrier[~ppl_carrier.bus.isin(n.buses.index)].copy()

                if not missing.empty:
                    bus_coords = n.buses.loc[ds.indexes["bus"], ["x", "y"]]

                    reassigned = 0
                    for _, row in missing.iterrows():
                        if pd.isna(row.lon) or pd.isna(row.lat):
                            continue

                        dist = (bus_coords.x - row.lon) ** 2 + (
                            bus_coords.y - row.lat
                        ) ** 2
                        row.bus = dist.idxmin()
                        valid = pd.concat(
                            [valid, pd.DataFrame([row])], ignore_index=False
                        )
                        reassigned += 1

                    if reassigned > 0:
                        logger.info(
                            f"{carrier}: reassigned {reassigned} plants without valid bus."
                        )

                caps_existing = (
                    valid.groupby("bus")["p_nom"]
                    .sum()
                    .reindex(ds.indexes["bus"], fill_value=0.0)
                )
            else:
                caps_existing = pd.Series(0.0, index=ds.indexes["bus"])

            p_max_pu = ds["profile"].transpose("time", "bus").to_pandas()
            p_nom_max = ds["p_nom_max"].to_pandas()
            weight = ds["weight"].to_pandas()

            # Gap filling using IRENA
            caps_final = caps_existing.copy()

            targets = get_irena_targets_for_carrier(
                carrier=carrier,
                estimate_renewable_capacities_config=estimate_cfg,
                countries=countries,
            )

            if not targets.empty:
                buses_country = n.buses.loc[ds.indexes["bus"], "country"]
                existing_by_country = caps_existing.groupby(buses_country).sum()
                gap_by_country = targets.sub(existing_by_country, fill_value=0.0).clip(
                    lower=0.0
                )

                if gap_by_country.sum() > 0:
                    for country, gap in gap_by_country.items():
                        if gap <= 0:
                            continue

                        buses_c = buses_country[buses_country == country].index
                        if len(buses_c) == 0:
                            logger.warning(
                                f"{carrier}: no buses found for country {country}, "
                                f"cannot redistribute {gap/1e3:.2f} GW from IRENA."
                            )
                            continue

                        add_per_bus = gap / len(buses_c)
                        caps_final.loc[buses_c] += add_per_bus

                existing_total = caps_existing.sum()
                target_total = targets.sum()
                gap_total = caps_final.sum() - caps_existing.sum()

                logger.info(
                    f"{carrier}: existing from powerplants = {existing_total/1e3:.2f} GW, "
                    f"IRENA target = {target_total/1e3:.2f} GW, "
                    f"gap filled = {gap_total/1e3:.2f} GW "
                    f"({100 * gap_total / max(target_total, 1):.1f}% of total)."
                )

            # Add generators
            suffix = " " + carrier

            renewable_lifetime = (
                costs.at["offwind", "lifetime"]
                if supcarrier == "offwind"
                else costs.at[carrier, "lifetime"]
            )

            n.madd(
                "Generator",
                ds.indexes["bus"],
                suffix,
                bus=ds.indexes["bus"],
                carrier=carrier,
                p_nom=caps_final,
                p_nom_min=caps_existing,
                p_nom_extendable=carrier in extendable_carriers["Generator"],
                p_nom_max=p_nom_max,
                p_max_pu=p_max_pu,
                weight=weight,
                marginal_cost=costs.at[supcarrier, "marginal_cost"],
                capital_cost=capital_cost,
                efficiency=costs.at[supcarrier, "efficiency"],
                lifetime=renewable_lifetime,
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
) -> None:
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

    logger.info(
        "Adding {} existing generators with capacities [GW] \n{}".format(
            len(ppl), ppl.groupby("carrier").p_nom.sum().div(1e3).round(2)
        )
    )

    n.madd(
        "Generator",
        ppl_grouped["bus"] + " " + ppl_grouped["carrier_gy"],
        carrier=ppl_grouped["carrier"],
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
    extendable_conventional = set(extendable_carriers["Generator"]) - set(
        renewable_carriers
    )

    if extendable_conventional:
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

        logger.info(f"Added extendable {extendable_conventional} generators.")
    else:
        logger.info("No extendable conventional generators configured.")

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


def apply_nuclear_p_max_pu(n, nuclear_p_max_pu):
    """
    Apply country-level static nuclear p_max_pu limits
    based on historical Energy Availability Factor (IAEA, 2022–2024).

    - If country is in CSV: apply p_max_pu
    - If country is NOT in CSV: keep default p_max_pu = 1.0 and warn
    """

    factors = nuclear_p_max_pu.set_index("country")["factor"].div(100.0)

    gens = n.generators[n.generators.carrier.str.startswith("nuclear", na=False)]
    if gens.empty:
        logger.info("No nuclear generators found: skipping nuclear p_max_pu limits.")
        return

    countries = gens.bus.map(n.buses.country)
    values = countries.map(factors)

    # Countries in model but missing in CSV
    missing = countries[values.isna()].dropna().unique()
    if len(missing) > 0:
        logger.warning(
            "No nuclear p_max_pu data for countries %s: "
            "keeping default p_max_pu = 1.0 for their nuclear generators.",
            ", ".join(sorted(missing)),
        )

    valid = values.notna()
    if not valid.any():
        logger.warning(
            "Nuclear p_max_pu CSV provided, but no matching countries found in the model."
        )
        return

    # Apply static constraint
    n.generators.loc[values.index[valid], "p_max_pu"] = values[valid]

    logger.info(
        "Applied static nuclear p_max_pu to %d nuclear generator(s) "
        "(source: IAEA 2022–2024).",
        valid.sum(),
    )


def attach_hydro(
    n: pypsa.Network,
    costs: pd.DataFrame,
    ppl: pd.DataFrame,
    hydro_min_inflow_pu: float = 1.0,
) -> None:
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
    tbd = ppl[ppl.technology.isna()]  # To be determined technologies

    inflow_idx = ror.index.union(hydro.index).union(tbd.index)
    if not inflow_idx.empty:
        with xr.open_dataarray(snakemake.input.profile_hydro) as inflow:
            found_plants = ppl.ppl_id[ppl.ppl_id.isin(inflow.indexes["plant"])]
            missing_plants_idxs = ppl.index.difference(found_plants.index)

            # if missing time series are found, notify the user and exclude missing hydro plants
            if not missing_plants_idxs.empty:
                # original total p_nom
                total_p_nom = ror.p_nom.sum() + hydro.p_nom.sum() + tbd.p_nom.sum()

                ror = ror.loc[ror.index.intersection(found_plants.index)]
                hydro = hydro.loc[hydro.index.intersection(found_plants.index)]
                tbd = tbd.loc[tbd.index.intersection(found_plants.index)]

                # loss of p_nom
                loss_p_nom = (
                    ror.p_nom.sum() + hydro.p_nom.sum() + tbd.p_nom.sum() - total_p_nom
                )

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
                inflow_agg = aggregate_inflow_by_group(ppl, ppl_grouped, inflow_t)

    # Heuristics for missing hydro technologies
    if not tbd.empty:
        inflow_pu_limit = inflow_t[tbd.index].mean() / tbd["p_nom"]  # Average MWh/MW
        mask_reservoir = inflow_pu_limit >= hydro_min_inflow_pu
        mask_ror = ~mask_reservoir
        to_be_hydro = tbd[mask_reservoir].copy()
        to_be_ror = tbd[mask_ror].copy()

        # Aggregate to_be_ror and to_be_hydro by (bus, carrier, grouping_year)
        to_be_hydro.loc[:, "carrier"] = "hydro"
        to_be_ror.loc[:, "carrier"] = "ror"
        to_be_ror_grouped = aggregate_ppl_by_bus_carrier_year(to_be_ror)
        to_be_hydro_grouped = aggregate_ppl_by_bus_carrier_year(to_be_hydro)
        inflow_agg_ror = aggregate_inflow_by_group(
            to_be_ror, to_be_ror_grouped, inflow_t
        )
        inflow_agg_hydro = aggregate_inflow_by_group(
            to_be_hydro, to_be_hydro_grouped, inflow_t
        )

        # Concatenate to existing ror and hydro dataframes and re-aggregate
        # to correctly merge any duplicate (bus, carrier_gy) groups
        ror = aggregate_ppl_by_bus_carrier_year(
            pd.concat(
                [
                    ppl[ppl["carrier"] == "ror"],  # original ungrouped ror plants
                    to_be_ror,  # ungrouped tbd reclassified as ror
                ]
            )
        )
        hydro = aggregate_ppl_by_bus_carrier_year(
            pd.concat(
                [
                    ppl[ppl["carrier"] == "hydro"],  # original ungrouped hydro plants
                    to_be_hydro,  # ungrouped tbd reclassified as hydro
                ]
            )
        )

        # Merge inflow: sum overlapping columns, add new ones
        inflow_agg = (
            pd.concat([inflow_agg, inflow_agg_ror, inflow_agg_hydro], axis=1)
            .groupby(level=0, axis=1)
            .sum()
        )

        logger.info(
            f"Identified {len(tbd)} hydro powerplants with unknown technology.\n"
            f"Hydropower plants with energy-to-capacity ratio ≥ {hydro_min_inflow_pu} "
            f"are classified as 'Reservoir'. The rest are 'Run-Of-River'.\n"
            f"Reservoir: {mask_reservoir.sum()} \n"
            f"Run-Of-River: {mask_ror.sum()}"
        )

    if "ror" in carriers and not ror.empty:
        n.madd(
            "Generator",
            ror.index,
            carrier=ror["carrier"],
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
            carrier=phs["carrier"],
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
        hydro_max_hours_default = c.get("hydro_max_hours_default", 6.0)
        if not missing_countries.empty:
            logger.warning(
                f"Assuming max_hours={hydro_max_hours_default} for hydro reservoirs in the countries: "
                "{}".format(", ".join(missing_countries))
            )
        hydro_max_hours = hydro.max_hours.where(
            hydro.max_hours > 0, hydro.country.map(max_hours_country)
        ).fillna(hydro_max_hours_default)

        n.madd(
            "StorageUnit",
            hydro.index,
            carrier=hydro["carrier"],
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


def attach_existing_batteries(
    n: pypsa.Network,
    costs: pd.DataFrame,
    ppl: pd.DataFrame,
    extendable_carriers: dict,
    max_hours: float,
) -> None:
    """
    Add existing battery storage units from powerplants.csv to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to modify.
    costs : pd.DataFrame
        DataFrame containing technology costs.
    ppl : pd.DataFrame
        Power plant DataFrame.
    extendable_carriers : dict
        Dictionary of extendable carriers for different component types.
    max_hours: float
        Amount of time it takes to fully charge batteries from empty if done at maximum power rate.

    Returns
    -------
    None
    """
    batteries = ppl.query('carrier == "battery"')
    if batteries.empty:
        logger.info("No existing batteries found in powerplants.csv.")
        return

    _add_missing_carriers_from_costs(n, costs, ["battery"])

    # Aggregate batteries by (bus, carrier, grouping_year)
    batteries_grouped = aggregate_ppl_by_bus_carrier_year(batteries)

    if "battery" in extendable_carriers["Store"]:
        battery_def = "Stores and Links"

        n.madd(
            "Bus",
            batteries_grouped["bus"] + " battery",
            location=batteries_grouped["bus"],
            carrier="battery",
            x=n.buses.loc[list(batteries_grouped["bus"])].x.values,
            y=n.buses.loc[list(batteries_grouped["bus"])].y.values,
        )

        # Add Stores for existing battery energy capacity
        n.madd(
            "Store",
            batteries_grouped.index,
            bus=batteries_grouped["bus"] + " battery",
            e_cyclic=True,
            e_nom_extendable=False,
            e_nom=batteries_grouped["p_nom"] * max_hours,
            carrier="battery",
            capital_cost=costs.at["battery storage", "capital_cost"],
            build_year=batteries_grouped["build_year"],
            lifetime=batteries_grouped["lifetime"],
        )

        # Add charger Links for existing battery power capacity
        n.madd(
            "Link",
            batteries_grouped.index.map(
                lambda x: x.replace(" battery", " battery charger")
            ),
            bus0=batteries_grouped["bus"].values,
            bus1=(batteries_grouped["bus"] + " battery").values,
            p_nom=batteries_grouped["p_nom"].values,
            p_nom_extendable=False,
            carrier="battery charger",
            efficiency=costs.at["battery inverter", "efficiency"] ** 0.5,
            capital_cost=costs.at["battery inverter", "capital_cost"],
            lifetime=costs.at["battery inverter", "lifetime"],
            build_year=batteries_grouped["build_year"].values,
        )

        # Add discharger Links for existing battery power capacity
        n.madd(
            "Link",
            batteries_grouped.index.map(
                lambda x: x.replace(" battery", " battery discharger")
            ),
            bus0=(batteries_grouped["bus"] + " battery").values,
            bus1=batteries_grouped["bus"].values,
            p_nom=batteries_grouped["p_nom"].values,
            p_nom_extendable=False,
            carrier="battery discharger",
            efficiency=costs.at["battery inverter", "efficiency"] ** 0.5,
            marginal_cost=costs.at["battery", "marginal_cost"],
            lifetime=costs.at["battery inverter", "lifetime"],
            build_year=batteries_grouped["build_year"].values,
        )

    else:
        battery_def = "StorageUnit"

        n.madd(
            "StorageUnit",
            batteries_grouped.index,
            bus=batteries_grouped["bus"],
            carrier=batteries_grouped["carrier"],
            p_nom=batteries_grouped["p_nom"],
            p_nom_extendable=False,
            p_nom_min=batteries_grouped["p_nom"],
            p_nom_max=batteries_grouped["p_nom"],
            capital_cost=costs.at["battery", "capital_cost"],
            max_hours=max_hours,
            efficiency_store=np.sqrt(costs.at["battery", "efficiency"]),
            efficiency_dispatch=np.sqrt(costs.at["battery", "efficiency"]),
            cyclic_state_of_charge=True,
            marginal_cost=costs.at["battery", "marginal_cost"],
            build_year=batteries_grouped["build_year"],
            lifetime=batteries_grouped["lifetime"],
        )

    logger.info(
        f"Added {len(batteries_grouped)} existing batteries defined as {battery_def} with total capacity "
        f"{batteries_grouped.p_nom.sum()/1e3:.2f} GW (max_hours={max_hours})."
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


def add_nice_carrier_names(n: pypsa.Network, config: dict) -> None:
    """
    Add nice names and colors to carriers.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to modify.
    config : dict
        The configuration dictionary containing plotting settings.

    Returns
    -------
    None
    """
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
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("add_electricity")

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.base_network)
    Nyears = n.snapshot_weightings.objective.sum() / 8760.0

    # Snakemake imports:
    demand_profiles = snakemake.input["demand_profiles"]

    costs = pd.read_csv(snakemake.input.tech_costs, index_col=0)
    ppl = load_powerplants(
        snakemake.input.powerplants,
        costs,
        snakemake.params.fill_values,
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
    attach_hydro(
        n, costs, ppl, snakemake.params.renewable["hydro"]["hydro_min_inflow_pu"]
    )

    attach_existing_batteries(
        n,
        costs,
        ppl,
        extendable_carriers,
        snakemake.params.electricity["max_hours"]["battery"],
    )
    apply_nuclear_p_max_pu(
        n,
        pd.read_csv(snakemake.input.nuclear_p_max_pu),
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

    # Log total installed capacities by carrier (GW)
    gen_caps = (
        n.generators.groupby(n.generators.carrier.str.strip().str.lower())
        .p_nom.sum()
        .div(1e3)
        .rename("Generators [GW]")
    )
    sto_caps = (
        n.storage_units.groupby(n.storage_units.carrier.str.strip().str.lower())
        .p_nom.sum()
        .div(1e3)
        .rename("StorageUnits [GW]")
    )

    summary = pd.concat([gen_caps, sto_caps], axis=1).fillna(0).sort_index()

    if not summary.empty:
        logger.info("\nInstalled capacities summary\n%s", summary.round(2))
    else:
        logger.info("No generators or storage units found.")

    n.export_to_netcdf(snakemake.output[0])
