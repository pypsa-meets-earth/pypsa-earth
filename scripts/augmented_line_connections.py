# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2021 The PyPSA-Africa Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later
# coding: utf-8
"""
Guarantees that every bus has at least X-number of connections.

Relevant Settings
-----------------

.. code:: yaml


.. seealso::


Inputs
------


Outputs
-------



Description
-----------


"""
import logging
import os

import networkx as nx
import numpy as np
import pandas as pd
import pypsa
from _helpers import configure_logging
from add_electricity import load_costs
from base_network import _set_links_underwater_fraction
from networkx.algorithms import complement
from networkx.algorithms.connectivity.edge_augmentation import k_edge_augmentation
from pypsa.geo import haversine_pts

logger = logging.getLogger(__name__)


# Functions
def haversine(p):
    coord0 = n.buses.loc[p.bus0, ["x", "y"]].values
    coord1 = n.buses.loc[p.bus1, ["x", "y"]].values
    return 1.5 * haversine_pts(coord0, coord1)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "augmented_line_connections", network="elec", simpl="", clusters="10"
        )
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)
    Nyears = n.snapshot_weightings.sum().values[0] / 8760.0
    costs = load_costs(
        Nyears,
        snakemake.input.tech_costs,
        snakemake.config["costs"],
        snakemake.config["electricity"],
    )
    options = snakemake.config["augmented_line_connection"]
    min_expansion_option = options.get("min_expansion")
    k_edge_option = options.get("connectivity_upgrade", 3)
    line_type_option = options.get("new_line_type", "HVDC")

    # k_edge algorithm implementation
    G = nx.Graph()

    network_buses = n.buses.loc[n.buses.carrier == "AC"].index
    G.add_nodes_from(np.unique(network_buses.values))

    # TODO: Currently only AC lines are read in and meshed. One need to combine
    # AC & DC lines and then move on.
    network_lines = n.lines
    sel = network_lines.s_nom > 100  # TODO: Check, should be all selected or filtered?
    attrs = ["bus0", "bus1", "length"]
    G.add_weighted_edges_from(network_lines.loc[:, attrs].values)

    # find all complement edges
    complement_edges = pd.DataFrame(complement(G).edges, columns=["bus0", "bus1"])
    complement_edges["length"] = complement_edges.apply(haversine, axis=1)

    # apply k_edge_augmentation weighted by length of complement edges
    k_edge = k_edge_option
    augmentation = k_edge_augmentation(G, k_edge, avail=complement_edges.values)
    new_network_lines = pd.DataFrame(augmentation, columns=["bus0", "bus1"])
    new_network_lines["length"] = new_network_lines.apply(haversine, axis=1)
    new_network_lines.index = new_network_lines.apply(
        lambda x: f"lines new {x.bus0} <-> {x.bus1}", axis=1
    )

    #  add new lines to the network
    if line_type_option == "HVDC":
        n.madd(
            "Link",
            new_network_lines.index,
            bus0=new_network_lines.bus0,
            bus1=new_network_lines.bus1,
            p_min_pu=-1,  # network is bidirectional
            p_nom_extendable=True,
            p_nom_min=min_expansion_option,
            length=new_network_lines.length,
            capital_cost=new_network_lines.length
            * costs.at["HVDC overhead", "capital_cost"],
            carrier="DC",
            lifetime=costs.at["HVDC overhead", "lifetime"],
            underwater_fraction=0.0,
        )

    if line_type_option == "HVAC":
        n.madd(
            "Line",
            new_network_lines.index,
            bus0=new_network_lines.bus0,
            bus1=new_network_lines.bus1,
            s_nom_extendable=True,
            # TODO: Check if minimum value needs to be set.
            s_nom_min=min_expansion_option,
            length=new_network_lines.length,
            capital_cost=new_network_lines.length
            * costs.at["HVAC overhead", "capital_cost"],
            carrier="AC",
            lifetime=costs.at["HVAC overhead", "lifetime"],
        )

        _set_links_underwater_fraction(n)

n.export_to_netcdf(snakemake.output.network)
