# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import pathlib
import sys

import pandas as pd
import pypsa

sys.path.append("./scripts")

from add_electricity import load_costs
from add_extra_components import (
    attach_hydrogen_pipelines,
    attach_storageunits,
    attach_stores,
)

path_cwd = pathlib.Path.cwd()
path_costs = pathlib.Path(path_cwd, "data", "costs.csv")

costs_config_dict = {
    "year": 2030,
    "version": "v0.5.0",
    "rooftop_share": 0.14,
    "USD2013_to_EUR2013": 0.7532,
    "fill_values": {
        "FOM": 0,
        "VOM": 0,
        "efficiency": 1,
        "fuel": 0,
        "investment": 0,
        "lifetime": 25,
        "CO2 intensity": 0,
        "discount rate": 0.07,
    },
    "marginal_cost": {
        "solar": 0.01,
        "onwind": 0.015,
        "offwind": 0.015,
        "hydro": 0.0,
        "H2": 0.0,
        "electrolysis": 0.0,
        "fuel cell": 0.0,
        "battery": 0.0,
        "battery inverter": 0.0,
    },
    "emission_prices": {
        "co2": 0.0,
    },
}

costs_config_df = pd.DataFrame.from_dict(costs_config_dict)

config_dict = {
    "electricity": {
        "base_voltage": 380.0,
        "voltages": [132.0, 220.0, 300.0, 380.0, 500.0, 750.0],
        "co2limit": 1.487e9,
        "co2base": 1.487e9,
        "agg_p_nom_limits": "data/agg_p_nom_minmax.csv",
        "hvdc_as_lines": True,
        "automatic_emission": True,
        "automatic_emission_base_year": 1990,
        "operational_reserve": {
            "activate": True,
            "epsilon_load": 0.02,
            "epsilon_vres": 0.02,
            "contingency": 0,
        },
        "max_hours": {
            "battery": 6,
            "H2": 168,
        },
        "extendable_carriers": {
            "Generator": ["solar", "onwind", "offwind-ac", "offwind-dc", "OCGT"],
            "StorageUnit": ["H2"],
            "Store": ["battery", "H2"],
            "Link": ["H2 pipeline"],
        },
        "powerplants_filter": "(DateOut >= 2022 or DateOut != DateOut)",
        "custom_powerplants": False,
        "conventional_carriers": [
            "nuclear",
            "oil",
            "OCGT",
            "CCGT",
            "coal",
            "lignite",
            "geothermal",
            "biomass",
        ],
        "renewable_carriers": [
            "solar",
            "csp",
            "onwind",
            "offwind-ac",
            "offwind-dc",
            "hydro",
        ],
        "estimate_renewable_capacities": {
            "stats": "irena",
            "year": 2020,
            "p_nom_min": 1,
            "p_nom_max": False,
            "technology_mapping": {
                "Offshore": ["offwind-ac", "offwind-dc"],
                "Onshore": ["onwind"],
                "PV": ["solar"],
            },
        },
    },
    "renewable": {
        "csp": {
            "cutout": "cutout-2013-era5-tutorial",
            "resource": {"method": "csp", "installation": "SAM_solar_tower"},
            "capacity_per_sqkm": 2.392,
            "copernicus": {
                "grid_codes": [20, 30, 40, 60, 90],
                "distancing_codes": [50],
                "distance_to_codes": 3000,
            },
            "natura": True,
            "potential": "simple",
            "clip_p_max_pu": 1.0e-2,
            "extendable": True,
            "csp_model": "simple",
        },
    },
}


def test_attach_storageunits():
    """
    Verify what is returned by attach_storageunits()
    """
    test_network_de = pypsa.examples.scigrid_de(from_master=True)
    number_years = test_network_de.snapshot_weightings.objective.sum() / 8760.0
    test_costs = load_costs(
        path_costs,
        costs_config_dict,
        config_dict["electricity"],
        number_years,
    )

    reference_component_dict = {
        "Bus": 585,
        "Carrier": 1,
        "Line": 852,
        "LineType": 34,
        "Transformer": 96,
        "TransformerType": 14,
        "Load": 489,
        "Generator": 1423,
        "StorageUnit": 623,
    }
    attach_storageunits(test_network_de, test_costs, config_dict)

    output_component_dict = {}
    for c in test_network_de.iterate_components(
        list(test_network_de.components.keys())[2:]
    ):
        output_component_dict[c.name] = len(c.df)

    assert output_component_dict == reference_component_dict


def test_attach_stores():
    """
    Verify what is returned by attach_stores()
    """
    test_network_de = pypsa.examples.scigrid_de(from_master=True)
    number_years = test_network_de.snapshot_weightings.objective.sum() / 8760.0
    test_costs = load_costs(
        path_costs,
        costs_config_dict,
        config_dict["electricity"],
        number_years,
    )

    reference_component_dict = {
        "Bus": 1755,
        "Carrier": 2,
        "Line": 852,
        "LineType": 34,
        "Transformer": 96,
        "TransformerType": 14,
        "Link": 2340,
        "Load": 489,
        "Generator": 1423,
        "StorageUnit": 38,
        "Store": 1170,
    }
    test_network_de.buses["country"] = "DE"
    attach_stores(test_network_de, test_costs, config_dict)

    output_component_dict = {}
    for c in test_network_de.iterate_components(
        list(test_network_de.components.keys())[2:]
    ):
        output_component_dict[c.name] = len(c.df)

    assert output_component_dict == reference_component_dict


def test_attach_hydrogen_pipelines():
    """
    Verify what is returned by attach_hydrogen_pipelines()
    """
    test_network_de = pypsa.examples.scigrid_de(from_master=True)
    number_years = test_network_de.snapshot_weightings.objective.sum() / 8760.0
    test_costs = load_costs(
        path_costs,
        costs_config_dict,
        config_dict["electricity"],
        number_years,
    )

    reference_component_dict = {
        "Bus": 585,
        "Line": 852,
        "LineType": 34,
        "Transformer": 96,
        "TransformerType": 14,
        "Link": 705,
        "Load": 489,
        "Generator": 1423,
        "StorageUnit": 38,
    }
    attach_hydrogen_pipelines(test_network_de, test_costs, config_dict)

    output_component_dict = {}
    for c in test_network_de.iterate_components(
        list(test_network_de.components.keys())[2:]
    ):
        output_component_dict[c.name] = len(c.df)

    assert output_component_dict == reference_component_dict
