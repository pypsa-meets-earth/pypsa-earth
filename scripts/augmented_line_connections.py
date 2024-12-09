# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
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
import os

import networkx as nx
import numpy as np
import pandas as pd
import pypsa
from _helpers import configure_logging, create_logger
from add_electricity import load_costs
from networkx.algorithms import complement
from networkx.algorithms.connectivity.edge_augmentation import k_edge_augmentation
from pypsa.geo import haversine_pts

logger = create_logger(__name__)


# Functions
def haversine(p):
    coord0 = n.buses.loc[p.bus0, ["x", "y"]].values
    coord1 = n.buses.loc[p.bus1, ["x", "y"]].values
    return 1.5 * haversine_pts(coord0, coord1)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "augmented_line_connections", network="elec", simpl="", clusters="4"
        )

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)
    Nyears = n.snapshot_weightings.sum().values[0] / 8760.0
    costs = load_costs(
        snakemake.input.tech_costs,
        snakemake.params.costs,
        snakemake.params.electricity,
        Nyears,
    )
    # TODO: Implement below comment in future. Requires transformer consideration.
    # component_type = {"False": "Line", "True":  "Link"}.get(snakemake.params.hvdc_as_lines)
    options = snakemake.params.augmented_line_connection
    min_expansion_option = options.get("min_expansion")
    k_edge_option = options.get("connectivity_upgrade", 3)
    line_type_option = options.get("new_line_type", ["HVDC"])
    min_DC_length = options.get("min_DC_length")

    # k_edge algorithm implementation
    G = nx.Graph()

    network_buses = n.buses.loc[n.buses.carrier == "AC"].index
    G.add_nodes_from(np.unique(network_buses.values))

    # TODO: Currently only AC lines are read in and meshed. One need to combine
    # AC & DC lines and then move on.
    # rather there is a need to filter-out DC lines
    network_lines = n.lines
    sel = network_lines.s_nom > 100  # TODO: Check, should be all selected or filtered?
    attrs = ["bus0", "bus1", "length"]
    G.add_weighted_edges_from(network_lines.loc[:, attrs].values)

    # find all complement edges, info on complement edges https://www.geeksforgeeks.org/complement-of-graph/
    complement_edges = pd.DataFrame(complement(G).edges, columns=["bus0", "bus1"])
    complement_edges["length"] = complement_edges.apply(haversine, axis=1)
    complement_edges["interconnector"] = np.invert(
        [
            (
                complement_edges.loc[x, "bus0"][0:2]
                == complement_edges.loc[x, "bus1"][0:2]
            )
            for x in complement_edges.index
        ]
    )

    # apply k_edge_augmentation weighted by length of complement edges
    # pick shortest lines per bus to fill k_edge condition (=degree of connectivity)
    k_edge = k_edge_option
    augmentation = k_edge_augmentation(
        G, k_edge, avail=complement_edges[["bus0", "bus1", "length"]].values
    )
    new_kedge_lines = pd.DataFrame(augmentation, columns=["bus0", "bus1"])
    new_kedge_lines["length"] = new_kedge_lines.apply(haversine, axis=1)
    new_kedge_lines.index = new_kedge_lines.apply(
        lambda x: f"lines new {x.bus0} <-> {x.bus1}", axis=1
    )

    # random sampling for long lines above <min DC length [km]>, including min and max distance, excluding interconnectors
    intracountry_edges = complement_edges[~complement_edges["interconnector"]]
    df = intracountry_edges[intracountry_edges["length"] > min_DC_length].sort_values(
        by=["length"]
    )
    random_sample = df.sample(
        frac=0.01, random_state=1
    )  # frac extract 1% of complement_edges as samples
    min_sample = df.head(1)
    max_sample = df.tail(1)
    new_long_lines = (
        pd.concat([min_sample, max_sample, random_sample]).drop_duplicates().dropna()
    )

    #  add new lines to the network
    if "HVDC" in list(line_type_option):
        n.madd(
            "Link",
            new_long_lines.index,
            suffix=" DC",
            bus0=new_long_lines.bus0,
            bus1=new_long_lines.bus1,
            type=snakemake.params.lines.get("dc_types"),
            p_min_pu=-1,  # network is bidirectional
            p_nom_extendable=True,
            p_nom_min=min_expansion_option,
            length=new_long_lines.length,
            capital_cost=new_long_lines.length
            * costs.at["HVDC overhead", "capital_cost"],
            carrier="DC",
            lifetime=costs.at["HVDC overhead", "lifetime"],
            underwater_fraction=0.0,
        )

    if "HVAC" in list(line_type_option):
        n.madd(
            "Line",
            new_kedge_lines.index,
            suffix=" AC",
            bus0=new_kedge_lines.bus0,
            bus1=new_kedge_lines.bus1,
            type=snakemake.params.lines["ac_types"].get(380),
            s_nom_extendable=True,
            # TODO: Check if minimum value needs to be set.
            s_nom_min=min_expansion_option,
            length=new_kedge_lines.length,
            capital_cost=new_kedge_lines.length
            * costs.at["HVAC overhead", "capital_cost"],
            carrier="AC",
            lifetime=costs.at["HVAC overhead", "lifetime"],
        )

        # TODO There is a need to add calculations of `underwater_fraction`
        # considering only lines added during augmentation
        # _set_dc_underwater_fraction(n.links, snakemake.input.regions_offshore)
        # _set_dc_underwater_fraction(n.lines, snakemake.input.regions_offshore)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output.network)
