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

- ``networks/elec_s{simpl}_{clusters}_pre_augmentation.nc``: Input network before augmentation.
- ``resources/costs_{year}.csv``: Technology cost assumptions.
- ``data/line_type_mapping_ac.csv``: Country-specific AC voltage-to-line-type mappings.
- ``data/line_type_mapping_dc.csv``: Country-specific DC voltage-to-line-type mappings.


Outputs
-------

- ``networks/elec_s{simpl}_{clusters}.nc``: Network with the configured connectivity augmentation applied.


Description
-----------

New HVAC connections use the closest country-specific line type available for the country and nominal voltage of their first bus. The default mapping is used when no country-specific mapping is available.
"""

import networkx as nx
import numpy as np
import pandas as pd
import pypsa
from _helpers import configure_logging, create_logger
from networkx.algorithms import complement
from networkx.algorithms.connectivity.edge_augmentation import k_edge_augmentation
from pypsa.geo import haversine_pts

logger = create_logger(__name__)


# Functions
def _load_linetypes_from_csv(path):
    """
    Load voltage-to-line-type mappings from a CSV file.
    """
    linetypes = pd.read_csv(path, index_col=0)
    linetypes.index = linetypes.index.astype(float)

    if "default" not in linetypes.columns:
        raise ValueError(f"Missing 'default' column in line type mapping file: {path}")

    return linetypes


def _get_linetype_by_voltage_and_country(v_nom, country, linetypes):
    """
    Return the closest available line type for a voltage and country.
    """
    if country in linetypes.columns:
        mapping = linetypes[country].dropna()
    else:
        mapping = pd.Series(dtype=object)

    if mapping.empty:
        mapping = linetypes["default"].dropna()

    if mapping.empty:
        raise ValueError(
            f"No line type mapping found for voltage {v_nom} kV "
            f"and country '{country}'."
        )

    voltage = min(mapping.index, key=lambda value: abs(value - v_nom))
    return mapping.at[voltage]


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

    clusters = str(getattr(snakemake.wildcards, "clusters", ""))
    ac_buses = n.buses.index[n.buses.carrier == "AC"]
    if clusters == "1" or len(ac_buses) < 2:
        logger.info(
            f"Skipping augmentation (clusters={clusters}, AC buses={len(ac_buses)})."
        )
        n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
        n.export_to_netcdf(snakemake.output.network)
        raise SystemExit(0)

    Nyears = n.snapshot_weightings.sum().values[0] / 8760.0
    costs = pd.read_csv(snakemake.input.tech_costs, index_col=0)
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
    ac_linetypes = _load_linetypes_from_csv(snakemake.input.line_type_mapping_ac)
    dc_linetypes = _load_linetypes_from_csv(snakemake.input.line_type_mapping_dc)

    dc_linetype = _get_linetype_by_voltage_and_country(
        500,
        None,
        dc_linetypes,
    )

    new_kedge_lines["country"] = new_kedge_lines["bus0"].map(n.buses["country"])
    new_kedge_lines["v_nom"] = new_kedge_lines["bus0"].map(n.buses["v_nom"])
    new_kedge_lines["type"] = new_kedge_lines.apply(
        lambda line: _get_linetype_by_voltage_and_country(
            line.v_nom,
            line.country,
            ac_linetypes,
        ),
        axis=1,
    )

    if "HVDC" in list(line_type_option):
        n.madd(
            "Link",
            new_long_lines.index,
            suffix=" DC",
            bus0=new_long_lines.bus0,
            bus1=new_long_lines.bus1,
            type=dc_linetype,
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
            type=new_kedge_lines["type"],
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
