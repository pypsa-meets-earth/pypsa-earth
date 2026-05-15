# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import sys

import pypsa
import pytest

sys.path.append("./scripts")
from test.conftest import battery_network

import _helpers
from solve_network import add_battery_constraints, add_chp_constraints


def test_add_battery_constraints(battery_network):

    # Apply the function
    add_battery_constraints(battery_network)

    constraints = battery_network.model.constraints["Link-charger_ratio"]
    assert constraints is not None

    # Check that the constraints were added
    assert "Link-charger_ratio" in battery_network.model.constraints

    # Check the constraint expression
    lhs = battery_network.model.constraints["Link-charger_ratio"][0].lhs

    assert lhs.equals(
        battery_network.model["Link-p_nom"].loc["battery charger"]
        - 0.9 * battery_network.model["Link-p_nom"].loc["battery discharger"]
    )


def test_chp_constraints(battery_network):

    # Apply the function
    add_chp_constraints(battery_network)

    constraints = battery_network.model.constraints["Link-CHP_ratio"]
    assert constraints is not None

    # Check that the constraints were added
    assert "Link-CHP_ratio" in battery_network.model.constraints

    # Check the constraint expression
    lhs = battery_network.model.constraints["Link-CHP_ratio"][0].lhs

    assert lhs.equals(
        battery_network.model["Link-p_nom"].loc["urban central CHP heat"]
        - 0.9 * battery_network.model["Link-p_nom"].loc["urban central CHP electric"]
    )
