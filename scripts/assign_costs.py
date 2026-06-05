# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Assign horizon-specific technology costs to the structurally-fixed electricity
network.

The electricity build pipeline is split into two phases:

1. **Structure** (built once per scenario, shared across all horizons):
   ``add_electricity`` → ``simplify_network`` → ``cluster_network`` →
   ``add_extra_components`` → ``elec_s{simpl}_{clusters}_ec.nc``

2. **Costs** (assigned once per planning horizon, cheap):
   ``assign_costs`` → ``elec_s{simpl}_{clusters}_ec_{planning_horizons}.nc``

This script performs phase 2.  It loads the per-horizon cost file and calls
``update_electricity_costs``, which refreshes every cost-vintage-dependent
attribute (capital_cost, marginal_cost, efficiency, lifetime) on all
components while leaving the structural attributes (topology, installed
capacities, build years, plant-specific efficiencies of existing units)
untouched.

Outputs
-------
- ``networks/elec_s{simpl}_{clusters}_ec_{planning_horizons}.nc``:
  Structurally identical to the input but with costs from the selected
  planning horizon.
"""

import pandas as pd
import pypsa
from _helpers import configure_logging, create_logger
from add_electricity import update_electricity_costs

_logger = create_logger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "assign_costs",
            simpl="",
            clusters="4",
            planning_horizons=2030,
        )

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)

    costs = pd.read_csv(snakemake.input.tech_costs, index_col=0)

    update_electricity_costs(
        n,
        costs,
        renewable_carriers=snakemake.params.electricity["renewable_carriers"],
        length_factor=snakemake.params.length_factor,
        hydro_capital_cost=snakemake.params.hydro_capital_cost,
    )

    _logger.info(
        "Assigned costs from planning horizon %s to network %s.",
        snakemake.wildcards.planning_horizons,
        snakemake.input.network,
    )

    n.export_to_netcdf(snakemake.output.network)
