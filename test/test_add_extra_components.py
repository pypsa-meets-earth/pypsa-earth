# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import pathlib
import sys

sys.path.append("./scripts")

from test.conftest import get_config_dict, get_power_network_scigrid_de

from add_electricity import load_costs
from add_extra_components import (
    attach_hydrogen_pipelines,
    attach_storageunits,
    attach_stores,
)

path_cwd = pathlib.Path.cwd()
path_costs = pathlib.Path(path_cwd, "data", "costs.csv")


def test_attach_storageunits(get_config_dict, get_power_network_scigrid_de):
    """
    Verify what is returned by attach_storageunits()
    """
    config_dict = get_config_dict
    config_dict["electricity"]["extendable_carriers"]["StorageUnit"] = ["H2"]
    test_network_de = get_power_network_scigrid_de
    number_years = test_network_de.snapshot_weightings.objective.sum() / 8760.0
    test_costs = load_costs(
        path_costs,
        config_dict["costs"],
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


def test_attach_stores(get_config_dict, get_power_network_scigrid_de):
    """
    Verify what is returned by attach_stores()
    """
    config_dict = get_config_dict
    config_dict["renewable"]["csp"]["csp_model"] = "simple"
    test_network_de = get_power_network_scigrid_de
    number_years = test_network_de.snapshot_weightings.objective.sum() / 8760.0
    test_costs = load_costs(
        path_costs,
        config_dict["costs"],
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


def test_attach_hydrogen_pipelines(get_config_dict, get_power_network_scigrid_de):
    """
    Verify what is returned by attach_hydrogen_pipelines()
    """
    config_dict = get_config_dict
    config_dict["electricity"]["extendable_carriers"]["Link"] = ["H2 pipeline"]
    test_network_de = get_power_network_scigrid_de
    number_years = test_network_de.snapshot_weightings.objective.sum() / 8760.0
    test_costs = load_costs(
        path_costs,
        config_dict["costs"],
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
