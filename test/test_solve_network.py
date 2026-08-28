# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Unit tests for the extra-functionality constraints in `solve_network.py`."""

import pandas as pd
import pypsa
from solve_network import add_battery_constraints, add_chp_constraints


def _base_network():
    """Return a minimal one-bus, three-snapshot network with a generator and a load."""
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=3, freq="h"))
    n.add("Bus", "bus0")
    n.add("Generator", "gen0", bus="bus0", p_nom=100, marginal_cost=10)
    n.add("Load", "load0", bus="bus0", p_set=50)
    return n


# -- add_chp_constraints ----------------------------------------------------


def test_add_chp_constraints_returns_when_no_links():
    """add_chp_constraints should return without raising on a network with no links."""
    n = _base_network()
    assert n.links.empty
    assert add_chp_constraints(n) is None


# -- add_battery_constraints ------------------------------------------------


def test_add_battery_constraints_returns_when_no_extendable_links():
    """add_battery_constraints should return without adding a constraint when no link is extendable."""
    n = _base_network()
    n.add(
        "Link",
        "battery charger",
        bus0="bus0",
        bus1="bus0",
        p_nom=10,
        p_nom_extendable=False,
        efficiency=0.95,
    )
    n.add(
        "Link",
        "battery discharger",
        bus0="bus0",
        bus1="bus0",
        p_nom=10,
        p_nom_extendable=False,
        efficiency=0.9,
    )
    n.optimize.create_model()
    assert add_battery_constraints(n) is None
    assert "Link-charger_ratio" not in n.model.constraints


def test_add_battery_constraints_adds_charger_ratio_constraint():
    """add_battery_constraints should add the Link-charger_ratio constraint when both sides are extendable."""
    n = _base_network()
    n.add(
        "Link",
        "battery charger",
        bus0="bus0",
        bus1="bus0",
        p_nom_extendable=True,
        efficiency=0.95,
    )
    n.add(
        "Link",
        "battery discharger",
        bus0="bus0",
        bus1="bus0",
        p_nom_extendable=True,
        efficiency=0.9,
    )
    n.optimize.create_model()
    add_battery_constraints(n)

    assert "Link-charger_ratio" in n.model.constraints
    constraint = n.model.constraints["Link-charger_ratio"]
    assert constraint.shape == (1,)
