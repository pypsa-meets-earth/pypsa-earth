# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Adds extra extendable components to the clustered and simplified network.

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
        extendable_carriers:
            StorageUnit:
            Store:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at :ref:`costs_cf`,
    :ref:`electricity_cf`

Inputs
------

- ``resources/costs.csv``: The database of cost assumptions for all included technologies for specific years from various sources; e.g. discount rate, lifetime, investment (CAPEX), fixed operation and maintenance (FOM), variable operation and maintenance (VOM), fuel costs, efficiency, carbon-dioxide intensity.

Outputs
-------

- ``networks/elec_s{simpl}_{clusters}_ec.nc``:


Description
-----------

The rule :mod:`add_extra_components` attaches additional extendable components to the clustered and simplified network. These can be configured in the ``config.yaml`` at ``electricity: extendable_carriers:``. It processes ``networks/elec_s{simpl}_{clusters}.nc`` to build ``networks/elec_s{simpl}_{clusters}_ec.nc``, which in contrast to the former (depending on the configuration) contain with **zero** initial capacity

- ``StorageUnits`` of carrier 'H2' and/or 'battery'. If this option is chosen, every bus is given an extendable ``StorageUnit`` of the corresponding carrier. The energy and power capacities are linked through a parameter that specifies the energy capacity as maximum hours at full dispatch power and is configured in ``electricity: max_hours:``. This linkage leads to one investment variable per storage unit. The default ``max_hours`` lead to long-term hydrogen and short-term battery storage units.

- ``Stores`` of carrier 'H2' and/or 'battery' in combination with ``Links``. If this option is chosen, the script adds extra buses with corresponding carrier where energy ``Stores`` are attached and which are connected to the corresponding power buses via two links, one each for charging and discharging. This leads to three investment variables for the energy capacity, charging and discharging capacity of the storage unit.
"""
import os

import numpy as np
import pandas as pd
import pypsa
from _helpers import (
    STORE_LOOKUP,
    configure_logging,
    create_logger,
    load_costs,
    lossy_bidirectional_links,
    override_component_attrs,
    set_length_based_efficiency,
)
from add_electricity import (
    _add_missing_carriers_from_costs,
    add_nice_carrier_names,
)

idx = pd.IndexSlice

logger = create_logger(__name__)


def get_available_storage_carriers(n, carriers, include_H2=False):
    """
    Filter and register available storage carriers from a given list.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network object to which valid storage carriers will be added.

    carriers : list of str
        A list of carrier names to filter and potentially register as storage technologies.

    Returns
    -------
    list of str
        A list of storage carriers that are both implemented and available in the network.

    Notes
    -----
    Hydrogen-related carriers (e.g., "H2", "H2 tank", "H2 underground") are
    currently excluded from the returned list, even if implemented.
    """
    implemented = set(STORE_LOOKUP.keys())
    input_carriers = set(carriers)

    not_implemented = input_carriers - implemented
    if not_implemented:
        logger.warning(
            f"The following carriers are not implemented as storage technologies in PyPSA-Eur: {sorted(not_implemented)}"
        )

    # Define carriers to exclude (e.g. handled differently elsewhere)

    excluded = set() if include_H2 else {"H2", "H2 tank", "H2 underground"}
    available_carriers = sorted((input_carriers & implemented) - excluded)

    # Add any missing carriers to the network
    missing_carriers = [c for c in available_carriers if c not in n.carriers.index]
    if missing_carriers:
        n.madd("Carrier", missing_carriers)

    return available_carriers


def attach_stores(n, costs, buses_i, carriers, include_H2=False, marginal_cost_storage=0):
    """
    Add storage technologies to a PyPSA network in the form of stores and links.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to which the storage components will be added.

    costs : pandas.DataFrame
        A DataFrame containing cost and technical parameters for each storage components.

    carriers : list of str
        A list of storage carrier names to be added to the network.
        Only those defined in `STORE_LOOKUP` and available in the network are used.
    """
    available_carriers = get_available_storage_carriers(n, carriers, include_H2)
    for carrier in available_carriers:
        roundtrip_correction = 0.5 if carrier in ["battery","li-ion"] else 1

        lookup = STORE_LOOKUP[carrier]
        lookup_store = lookup["store"]
        if "bicharger" in lookup:
            lookup_charge = lookup_discharge = lookup["bicharger"]
        else:
            lookup_charge = lookup["charger"]
            lookup_discharge = lookup["discharger"]

        bus_names = buses_i + f" {carrier}"

        n.madd(
            "Bus",
            bus_names,
            location=buses_i,
            carrier=carrier,
            x=n.buses.loc[list(buses_i)].x.values,
            y=n.buses.loc[list(buses_i)].y.values,
        )

        n.madd(
            "Store",
            bus_names,
            bus=bus_names,
            e_cyclic=True,
            e_nom_extendable=True,
            carrier=carrier,
            capital_cost=costs.at[lookup_store, "capital_cost"],
            lifetime=costs.at[lookup_store, "lifetime"],
        )

        n.madd(
            "Link",
            bus_names,
            suffix=" charger",
            bus0=buses_i,
            bus1=bus_names,
            carrier=f"{carrier} charger",
            efficiency=costs.at[lookup_charge, "efficiency"] 
            ** roundtrip_correction,
            capital_cost=costs.at[lookup_charge, "capital_cost"],
            p_nom_extendable=True,
            lifetime=costs.at[lookup_charge, "lifetime"],
        )

        n.madd(
            "Link",
            bus_names,
            suffix=" discharger",
            bus0=bus_names,
            bus1=buses_i,
            carrier=f"{carrier} discharger",
            efficiency=costs.at[lookup_discharge, "efficiency"] 
            ** roundtrip_correction,
            marginal_cost=marginal_cost_storage,
            p_nom_extendable=True,
            lifetime=costs.at[lookup_discharge, "lifetime"],
        )

    logger.info("Add the following technologies as stores and links:\n - " + "\n - ".join(available_carriers))


def attach_storageunits(n, costs, buses_i, carriers, max_hours, include_H2=False, marginal_cost_storage=0):
    """
    Add storage technologies to a PyPSA network in the form of storage units.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to which the storage components will be added.

    costs : pandas.DataFrame
        A DataFrame containing cost and technical parameters for each storage components.

    carriers : list of str
        A list of storage carrier names to be added to the network.
        Only those defined in `STORE_LOOKUP` and available in the network are used.

    max_hours : dict
        Dictionary mapping each carrier to its maximum storage duration in hours.
    """
    available_carriers = get_available_storage_carriers(n, carriers, include_H2)
    for carrier in available_carriers:
        max_hour = max_hours.get(carrier)
        if max_hour is None:
            logger.warning(f"No max_hours defined for carrier '{carrier}'. Skipping.")
            continue

        roundtrip_correction = 0.5 if carrier in ["battery","li-ion"] else 1

        lookup = STORE_LOOKUP[carrier]
        if "bicharger" in lookup:
            lookup_charge = lookup_discharge = lookup["bicharger"]
        else:
            lookup_charge = lookup["charger"]
            lookup_discharge = lookup["discharger"]

        n.madd(
            "StorageUnit",
            buses_i + f" {carrier}",
            bus=buses_i,
            carrier=carrier,
            p_nom_extendable=True,
            capital_cost=costs.at[carrier, "capital_cost"],
            marginal_cost=marginal_cost_storage,
            efficiency_store=costs.at[lookup_charge, "efficiency"]
            ** roundtrip_correction,
            efficiency_dispatch=costs.at[lookup_discharge, "efficiency"]
            ** roundtrip_correction,
            max_hours=max_hour,
            cyclic_state_of_charge=True,
            lifetime=costs.at[carrier, "lifetime"],
        )

    logger.info("Add the following technologies as storage units:\n - " + "\n - ".join(available_carriers))


def attach_advance_csp(n, costs):
    # add separate buses for csp
    main_buses = n.generators.query("carrier == 'csp'").bus
    csp_buses_i = n.madd(
        "Bus",
        main_buses + " csp",
        carrier="csp",
        x=n.buses.loc[main_buses, "x"].values,
        y=n.buses.loc[main_buses, "y"].values,
        country=n.buses.loc[main_buses, "country"].values,
    )
    n.generators.loc[main_buses.index, "bus"] = csp_buses_i

    # add stores for csp
    n.madd(
        "Store",
        csp_buses_i,
        bus=csp_buses_i,
        carrier="csp",
        e_cyclic=True,
        e_nom_extendable=True,
        capital_cost=costs.at["csp-tower TES", "capital_cost"],
        marginal_cost=costs.at["csp-tower TES", "marginal_cost"],
    )

    # add links for csp
    n.madd(
        "Link",
        csp_buses_i,
        bus0=csp_buses_i,
        bus1=main_buses,
        carrier="csp",
        efficiency=costs.at["csp-tower", "efficiency"],
        capital_cost=costs.at["csp-tower", "capital_cost"],
        p_nom_extendable=True,
        marginal_cost=costs.at["csp-tower", "marginal_cost"],
    )


def attach_hydrogen_pipelines(n, costs, config, transmission_efficiency):
    elec_opts = config["electricity"]
    ext_carriers = elec_opts["extendable_carriers"]
    as_stores = ext_carriers.get("Store", [])

    if "H2 pipeline" not in ext_carriers.get("Link", []):
        return

    assert "H2" in as_stores, (
        "Attaching hydrogen pipelines requires hydrogen "
        "storage to be modelled as Store-Link-Bus combination. See "
        "`config.yaml` at `electricity: extendable_carriers: Store:`."
    )

    # determine bus pairs
    attrs = ["bus0", "bus1", "length"]
    candidates = pd.concat(
        [n.lines[attrs], n.links.query('carrier=="DC"')[attrs]]
    ).reset_index(drop=True)

    # remove bus pair duplicates regardless of order of bus0 and bus1
    h2_links = candidates[
        ~pd.DataFrame(np.sort(candidates[["bus0", "bus1"]])).duplicated()
    ]
    h2_links.index = h2_links.apply(lambda c: f"H2 pipeline {c.bus0}-{c.bus1}", axis=1)

    # add pipelines
    n.madd(
        "Link",
        h2_links.index,
        bus0=h2_links.bus0.values + " H2",
        bus1=h2_links.bus1.values + " H2",
        p_min_pu=-1,
        p_nom_extendable=True,
        length=h2_links.length.values,
        capital_cost=costs.at["H2 pipeline", "capital_cost"] * h2_links.length,
        carrier="H2 pipeline",
    )

    # split the pipeline into two unidirectional links to properly apply transmission losses in both directions.
    lossy_bidirectional_links(n, "H2 pipeline")

    # set the pipelines efficiency and the electricity required by the pipeline for compression
    set_length_based_efficiency(n, "H2 pipeline", " H2", transmission_efficiency)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("add_extra_components", simpl="", clusters=10)

    configure_logging(snakemake)

    overrides = override_component_attrs(snakemake.input.overrides)
    n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)
    Nyears = n.snapshot_weightings.objective.sum() / 8760.0
    transmission_efficiency = snakemake.params.transmission_efficiency
    config = snakemake.config

    costs = load_costs(
        snakemake.input.tech_costs,
        config["costs"],
        config["costs"]["output_currency"],
        config["costs"]["fill_values"],
        config["electricity"]["max_hours"],
        Nyears,
        config["costs"]["default_exchange_rate"],
        config["costs"]["future_exchange_rate_strategy"],
        config["costs"]["custom_future_exchange_rate"],
    )

    buses_i = n.buses.index
    attach_stores(
        n, 
        costs, 
        buses_i,
        config["electricity"]["extendable_carriers"]["Store"],
        include_H2=True,
        marginal_cost_storage=config["sector"]["marginal_cost_storage"],
    )

    attach_storageunits(
        n, 
        costs,
        buses_i,
        config["electricity"]["extendable_carriers"]["StorageUnit"], 
        config["electricity"]["max_hours"],
        include_H2=True,
        marginal_cost_storage=config["sector"]["marginal_cost_storage"],
    )

    if ("csp" in config["electricity"]["renewable_carriers"]) and (
        config["renewable"]["csp"]["csp_model"] == "advanced"
    ):
        attach_advance_csp(n, costs)

    attach_hydrogen_pipelines(n, costs, config, transmission_efficiency)

    add_nice_carrier_names(n, config=snakemake.config)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output[0])
