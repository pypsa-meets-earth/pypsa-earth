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
    configure_logging,
    create_logger,
    lossy_bidirectional_links,
    override_component_attrs,
    set_length_based_efficiency,
)
from add_electricity import (
    _add_missing_carriers_from_costs,
    add_nice_carrier_names,
    load_costs,
)

idx = pd.IndexSlice

logger = create_logger(__name__)


def attach_storageunits(n, costs, config):
    carriers = config["electricity"]["extendable_carriers"]["StorageUnit"]
    carriers_classic = [x for x in carriers if x == "H2" or x == "battery"]
    carriers_database = [x for x in carriers if x != "H2" and x != "battery"]
    max_hours = config["electricity"]["max_hours"]

    _add_missing_carriers_from_costs(n, costs, carriers)

    buses_i = n.buses.index

    lookup_store = {"H2": "electrolysis", "battery": "battery inverter"}
    lookup_dispatch = {"H2": "fuel cell", "battery": "battery inverter"}

    # TODO: add metadata to H2 and battery to only read in over database interface
    for carrier in carriers_classic:
        roundtrip_correction = 0.5 if carrier == "battery" else 1
        n.madd(
            "StorageUnit",
            buses_i,
            " " + carrier,
            bus=buses_i,
            carrier=carrier,
            p_nom_extendable=True,
            capital_cost=costs.at[carrier, "capital_cost"],
            marginal_cost=costs.at[carrier, "marginal_cost"],
            efficiency_store=costs.at[lookup_store[carrier], "efficiency"]
            ** roundtrip_correction,
            efficiency_dispatch=costs.at[lookup_dispatch[carrier], "efficiency"]
            ** roundtrip_correction,
            max_hours=max_hours[carrier],
            cyclic_state_of_charge=True,
        )

    # automated attach of storage_units from database
    for carrier in carriers_database:
        tech_type = costs.technology_type
        charger_or_bicharger_filter = (costs.carrier == carrier) & (
            (tech_type == "charger") | (tech_type == "bicharger")
        )
        discharger_or_bicharger_filter = (costs.carrier == carrier) & (
            (tech_type == "discharger") | (tech_type == "bicharger")
        )
        n.madd(
            "StorageUnit",
            buses_i,
            " " + carrier,
            bus=buses_i,
            carrier=carrier,
            p_nom_extendable=True,
            capital_cost=costs.at[carrier, "capital_cost"],
            marginal_cost=costs.at[carrier, "marginal_cost"],
            efficiency_store=float(
                costs.loc[charger_or_bicharger_filter, "efficiency"]
            ),
            efficiency_dispatch=float(
                costs.loc[discharger_or_bicharger_filter, "efficiency"]
            ),
            max_hours=max_hours[carrier],
            cyclic_state_of_charge=True,
        )


def attach_stores(n, costs, config):
    """
    Add storage units as store and links.

    Two types of storage exists:
    1. where every components can be independently sized e.g. hydrogen
    2. where charger and discharger are the same component e.g. Li-battery (inverter)
    """
    carriers = config["electricity"]["extendable_carriers"]["Store"]
    carriers_classic = [x for x in carriers if x == "H2" or x == "battery"]
    carriers_database = [x for x in carriers if x != "H2" and x != "battery"]
    _add_missing_carriers_from_costs(n, costs, carriers)

    buses_i = n.buses.index
    bus_sub_dict = {k: n.buses[k].values for k in ["x", "y", "country"]}

    # TODO: add metadata to H2 and battery to only read in over database interface
    if "H2" in carriers_classic:
        h2_buses_i = n.madd("Bus", buses_i + " H2", carrier="H2", **bus_sub_dict)

        n.madd(
            "Store",
            h2_buses_i,
            bus=h2_buses_i,
            carrier="H2",
            e_nom_extendable=True,
            e_cyclic=True,
            capital_cost=costs.at["hydrogen storage tank", "capital_cost"],
        )

        n.madd(
            "Link",
            h2_buses_i + " Electrolysis",
            bus0=buses_i,
            bus1=h2_buses_i,
            carrier="H2 electrolysis",
            p_nom_extendable=True,
            efficiency=costs.at["electrolysis", "efficiency"],
            capital_cost=costs.at["electrolysis", "capital_cost"],
            marginal_cost=costs.at["electrolysis", "marginal_cost"],
        )

        n.madd(
            "Link",
            h2_buses_i + " Fuel Cell",
            bus0=h2_buses_i,
            bus1=buses_i,
            carrier="H2 fuel cell",
            p_nom_extendable=True,
            efficiency=costs.at["fuel cell", "efficiency"],
            capital_cost=costs.at["fuel cell", "capital_cost"],
            marginal_cost=costs.at["fuel cell", "marginal_cost"],
        )

    if "battery" in carriers_classic:
        b_buses_i = n.madd(
            "Bus", buses_i + " battery", carrier="battery", **bus_sub_dict
        )

        n.madd(
            "Store",
            b_buses_i,
            bus=b_buses_i,
            carrier="battery",
            e_cyclic=True,
            e_nom_extendable=True,
            capital_cost=costs.at["battery storage", "capital_cost"],
            marginal_cost=costs.at["battery", "marginal_cost"],
        )

        n.madd(
            "Link",
            b_buses_i + " charger",
            bus0=buses_i,
            bus1=b_buses_i,
            carrier="battery charger",
            # NB: the efficiencies are "round trip efficiencies"
            efficiency=costs.at["battery inverter", "efficiency"],
            capital_cost=costs.at["battery inverter", "capital_cost"],
            p_nom_extendable=True,
            marginal_cost=costs.at["battery inverter", "marginal_cost"],
        )

        n.madd(
            "Link",
            b_buses_i + " discharger",
            bus0=b_buses_i,
            bus1=buses_i,
            carrier="battery discharger",
            # NB: the efficiencies are "round trip efficiencies"
            efficiency=costs.at["battery inverter", "efficiency"],
            p_nom_extendable=True,
            marginal_cost=costs.at["battery inverter", "marginal_cost"],
        )

    # automated attach of store-link storages from database
    for c in carriers_database:
        carrier_buses_i = n.madd("Bus", buses_i + f" {c}", carrier=c, **bus_sub_dict)
        tech_type = costs.technology_type
        carrier = costs.carrier
        # filters for pypsa store
        store_filter = (carrier == c) & (tech_type == "store")
        # filters for pypsa charger or bicharger
        charger_or_bicharger_filter = (carrier == c) & (
            (tech_type == "charger") | (tech_type == "bicharger")
        )
        # filters for pypsa discharger or bicharger
        discharger_or_bicharger_filter = (carrier == c) & (
            (tech_type == "discharger") | (tech_type == "bicharger")
        )
        n.madd(
            "Store",
            carrier_buses_i,
            bus=carrier_buses_i,
            carrier=c,
            e_nom_extendable=True,
            e_cyclic=True,
            capital_cost=float(costs.loc[store_filter, "capital_cost"]),
        )

        n.madd(
            "Link",
            carrier_buses_i + " charger",
            bus0=buses_i,
            bus1=carrier_buses_i,
            carrier=f"{c} charger",
            p_nom_extendable=True,
            efficiency=float(costs.loc[charger_or_bicharger_filter, "efficiency"]),
            capital_cost=float(costs.loc[charger_or_bicharger_filter, "capital_cost"]),
            marginal_cost=float(
                costs.loc[charger_or_bicharger_filter, "marginal_cost"]
            ),
        )

        # capital cost of discharger is None if bicharger since already added in previous link
        if (
            costs.loc[discharger_or_bicharger_filter].technology_type.item()
            == "bicharger"
        ):
            zero_or_value = 0
        else:
            zero_or_value = float(
                costs.loc[discharger_or_bicharger_filter, "capital_cost"]
            )

        n.madd(
            "Link",
            carrier_buses_i + " discharger",
            bus0=carrier_buses_i,
            bus1=buses_i,
            carrier=f"{c} discharger",
            p_nom_extendable=True,
            efficiency=float(costs.loc[discharger_or_bicharger_filter, "efficiency"]),
            capital_cost=zero_or_value,
            marginal_cost=float(
                costs.loc[discharger_or_bicharger_filter, "marginal_cost"]
            ),
        )

    if ("csp" in config["electricity"]["renewable_carriers"]) and (
        config["renewable"]["csp"]["csp_model"] == "advanced"
    ):
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
    ext_carriers = config["electricity"]["extendable_carriers"]
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
        config["electricity"],
        Nyears,
    )

    attach_storageunits(n, costs, config)
    attach_stores(n, costs, config)
    attach_hydrogen_pipelines(n, costs, config, transmission_efficiency)

    add_nice_carrier_names(n, config=config)

    n.meta = dict(config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output[0])
