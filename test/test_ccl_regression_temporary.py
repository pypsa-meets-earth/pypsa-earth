# SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Temporary CCL regression tests for PR #2110."""

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pypsa
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import solve_network


@pytest.mark.parametrize(
    ("include_existing", "expected_min", "expected_max"),
    [
        (False, 50.0, 100.0),
        (True, 30.0, 80.0),
    ],
)
def test_ccl_with_fixed_and_extendable_generators(
    monkeypatch, include_existing, expected_min, expected_max
):
    """Reproduce the mismatch between static generators and model variables."""
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2030-01-01", periods=2, freq="h"))

    n.add("Bus", "DE0", country="DE")
    n.add("Carrier", "solar")

    n.add(
        "Generator",
        "solar_fixed",
        bus="DE0",
        carrier="solar",
        p_nom=20.0,
    )

    n.add(
        "Generator",
        "solar_extendable",
        bus="DE0",
        carrier="solar",
        p_nom_extendable=True,
        p_nom_max=200.0,
        capital_cost=1.0,
    )

    n.add("Load", "demand", bus="DE0", p_set=50.0)

    n.optimize.create_model()

    capacity_variable = n.model["Generator-p_nom"]
    generator_dim = capacity_variable.dims[-1]
    model_generators = capacity_variable.indexes[generator_dim]

    # Original CCL selection includes fixed and extendable generators.
    original_generators = n.generators[n.generators.carrier.eq("solar")].index

    assert set(original_generators) == {
        "solar_fixed",
        "solar_extendable",
    }

    # Only the extendable generator has a nominal-capacity variable.
    assert list(model_generators) == ["solar_extendable"]
    assert not original_generators.equals(model_generators)

    limits = pd.DataFrame(
        {
            ("2030", "min"): [50.0],
            ("2030", "max"): [100.0],
        },
        index=pd.MultiIndex.from_tuples(
            [("DE", "solar")],
            names=["country", "carrier"],
        ),
    )

    monkeypatch.setattr(
        solve_network,
        "snakemake",
        SimpleNamespace(
            input=SimpleNamespace(agg_p_nom_minmax="unused.csv"),
            wildcards=SimpleNamespace(planning_horizons="2030"),
        ),
        raising=False,
    )

    monkeypatch.setattr(
        solve_network,
        "read_csv_nafix",
        lambda *args, **kwargs: limits,
    )

    config = {
        "electricity": {
            "agg_p_nom_limits": {
                "include_existing": include_existing,
            }
        }
    }

    solve_network.add_CCL_constraints(n, config)

    minimum = np.asarray(n.model.constraints["agg_p_nom_min"].rhs).ravel()

    maximum = np.asarray(n.model.constraints["agg_p_nom_max"].rhs).ravel()

    assert minimum.tolist() == pytest.approx([expected_min])
    assert maximum.tolist() == pytest.approx([expected_max])
