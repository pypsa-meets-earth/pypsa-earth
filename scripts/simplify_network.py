# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Lifts electrical transmission network to a single 380 kV voltage layer, removes
dead-ends of the network, and reduces multi-hop HVDC connections to a single
link.

Relevant Settings
-----------------

.. code:: yaml

    clustering:
        simplify:
        aggregation_strategies:

    costs:
        year:
        version:
        rooftop_share:
        USD2013_to_EUR2013:
        dicountrate:
        emission_prices:

    electricity:
        max_hours:

    lines:
        length_factor:

    links:
        p_max_pu:

    solving:
        solver:
            name:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`costs_cf`, :ref:`electricity_cf`, :ref:`renewable_cf`,
    :ref:`lines_cf`, :ref:`links_cf`, :ref:`solving_cf`

Inputs
------

- ``resources/costs.csv``: The database of cost assumptions for all included technologies for specific years from various sources; e.g. discount rate, lifetime, investment (CAPEX), fixed operation and maintenance (FOM), variable operation and maintenance (VOM), fuel costs, efficiency, carbon-dioxide intensity.
- ``resources/regions_onshore.geojson``: confer :ref:`busregions`
- ``resources/regions_offshore.geojson``: confer :ref:`busregions`
- ``networks/elec.nc``: confer :ref:`electricity`

Outputs
-------

- ``resources/regions_onshore_elec_s{simpl}.geojson``:

    .. image:: /img/regions_onshore_elec_s.png
            :width: 33 %

- ``resources/regions_offshore_elec_s{simpl}.geojson``:

    .. image:: /img/regions_offshore_elec_s  .png
            :width: 33 %

- ``resources/busmap_elec_s{simpl}.csv``: Mapping of buses from ``networks/elec.nc`` to ``networks/elec_s{simpl}.nc``;
- ``networks/elec_s{simpl}.nc``:

    .. image:: /img/elec_s.png
        :width: 33 %

Description
-----------

The rule :mod:`simplify_network` does up to four things:

1. Create an equivalent transmission network in which all voltage levels are mapped to the 380 kV level by the function ``simplify_network(...)``.

2. DC only sub-networks that are connected at only two buses to the AC network are reduced to a single representative link in the function ``simplify_links(...)``. The components attached to buses in between are moved to the nearest endpoint. The grid connection cost of offshore wind generators are added to the capital costs of the generator.

3. Stub lines and links, i.e. dead-ends of the network, are sequentially removed from the network in the function ``remove_stubs(...)``. Components are moved along.

4. Optionally, if an integer were provided for the wildcard ``{simpl}`` (e.g. ``networks/elec_s500.nc``), the network is clustered to this number of clusters with the routines from the ``cluster_network`` rule with the function ``cluster_network.cluster(...)``. This step is usually skipped!
"""
import logging
import os
import sys
from functools import reduce

import numpy as np
import pandas as pd
import pypsa
import scipy as sp
from _helpers import configure_logging, get_aggregation_strategies, update_p_nom_max
from add_electricity import load_costs
from cluster_network import cluster_regions, clustering_for_n_clusters
from pypsa.io import import_components_from_dataframe, import_series_from_dataframe
from pypsa.networkclustering import (
    aggregategenerators,
    aggregateoneport,
    busmap_by_stubs,
    get_clustering_from_busmap,
)
from scipy.sparse.csgraph import connected_components, dijkstra

sys.settrace

logger = logging.getLogger(__name__)


def simplify_network_to_base_voltage(n, linetype, base_voltage):
    """
    Fix all lines to a voltage level of base voltage level and remove all
    transformers.

    The function preserves the transmission capacity for each line while
    updating its voltage level, line type and number of parallel bundles
    (num_parallel). Transformers are removed and connected components
    are moved from their starting bus to their ending bus. The
    corresponding starting buses are removed as well.
    """
    logger.info(f"Mapping all network lines onto a single {int(base_voltage)}kV layer")

    n.lines["type"] = linetype
    n.lines["v_nom"] = base_voltage
    n.lines["i_nom"] = n.line_types.i_nom[linetype]
    # Note: s_nom is set in base_network
    n.lines["num_parallel"] = n.lines.eval("s_nom / (sqrt(3) * v_nom * i_nom)")

    # Replace transformers by lines
    trafo_map = pd.Series(n.transformers.bus1.values, n.transformers.bus0.values)
    trafo_map = trafo_map[~trafo_map.index.duplicated(keep="first")]
    several_trafo_b = trafo_map.isin(trafo_map.index)
    trafo_map[several_trafo_b] = trafo_map[several_trafo_b].map(trafo_map)
    missing_buses_i = n.buses.index.difference(trafo_map.index)
    trafo_map = pd.concat([trafo_map, pd.Series(missing_buses_i, missing_buses_i)])

    for c in n.one_port_components | n.branch_components:
        df = n.df(c)
        for col in df.columns:
            if col.startswith("bus"):
                df[col] = df[col].map(trafo_map)

    n.mremove("Transformer", n.transformers.index)
    n.mremove("Bus", n.buses.index.difference(trafo_map))

    return n, trafo_map


def _prepare_connection_costs_per_link(n, costs, config):
    if n.links.empty:
        return {}

    connection_costs_per_link = {}

    for tech in config["renewable"]:
        if tech.startswith("offwind"):
            connection_costs_per_link[tech] = (
                n.links.length
                * config["lines"]["length_factor"]
                * (
                    n.links.underwater_fraction
                    * costs.at[tech + "-connection-submarine", "capital_cost"]
                    + (1.0 - n.links.underwater_fraction)
                    * costs.at[tech + "-connection-underground", "capital_cost"]
                )
            )

    return connection_costs_per_link


def _compute_connection_costs_to_bus(
    n, busmap, costs, config, connection_costs_per_link=None, buses=None
):
    if connection_costs_per_link is None:
        connection_costs_per_link = _prepare_connection_costs_per_link(n, costs, config)

    if buses is None:
        buses = busmap.index[busmap.index != busmap.values]

    connection_costs_to_bus = pd.DataFrame(index=buses)

    for tech in connection_costs_per_link:
        adj = n.adjacency_matrix(
            weights=pd.concat(
                dict(
                    Link=connection_costs_per_link[tech]
                    .reindex(n.links.index)
                    .astype(float),
                    Line=pd.Series(0.0, n.lines.index),
                )
            )
        )
        costs_between_buses = dijkstra(
            adj, directed=False, indices=n.buses.index.get_indexer(buses)
        )
        connection_costs_to_bus[tech] = costs_between_buses[
            np.arange(len(buses)), n.buses.index.get_indexer(busmap.loc[buses])
        ]
    return connection_costs_to_bus


def _adjust_capital_costs_using_connection_costs(n, connection_costs_to_bus, output):
    connection_costs = {}
    for tech in connection_costs_to_bus:
        tech_b = n.generators.carrier == tech
        costs = (
            n.generators.loc[tech_b, "bus"]
            .map(connection_costs_to_bus[tech])
            .loc[lambda s: s > 0]
        )
        if not costs.empty:
            n.generators.loc[costs.index, "capital_cost"] += costs
            logger.info(
                "Displacing {} generator(s) and adding connection costs to capital_costs: {} ".format(
                    tech,
                    ", ".join(
                        "{:.0f} Eur/MW/a for `{}`".format(d, b)
                        for b, d in costs.items()
                    ),
                )
            )
            connection_costs[tech] = costs
    pd.DataFrame(connection_costs).to_csv(output.connection_costs)


def _aggregate_and_move_components(
    n,
    busmap,
    connection_costs_to_bus,
    output,
    aggregate_one_ports={"Load", "StorageUnit"},
    aggregation_strategies=dict(),
    exclude_carriers=None,
):
    def replace_components(n, c, df, pnl):
        n.mremove(c, n.df(c).index)

        import_components_from_dataframe(n, df, c)
        for attr, df in pnl.items():
            if not df.empty:
                import_series_from_dataframe(n, df, c, attr)

    _adjust_capital_costs_using_connection_costs(n, connection_costs_to_bus, output)

    _, generator_strategies = get_aggregation_strategies(aggregation_strategies)

    carriers = set(n.generators.carrier) - set(exclude_carriers)
    generators, generators_pnl = aggregategenerators(
        n, busmap, carriers=carriers, custom_strategies=generator_strategies
    )

    replace_components(n, "Generator", generators, generators_pnl)

    for one_port in aggregate_one_ports:
        df, pnl = aggregateoneport(n, busmap, component=one_port)
        replace_components(n, one_port, df, pnl)

    buses_to_del = n.buses.index.difference(busmap)
    n.mremove("Bus", buses_to_del)
    for c in n.branch_components:
        df = n.df(c)
        n.mremove(c, df.index[df.bus0.isin(buses_to_del) | df.bus1.isin(buses_to_del)])


# Filter AC lines to avoid mixing with DC part when processing links
def contains_ac(ls):
    def is_ac_branch(x):
        if x[0] == "Link":
            return n.links.loc[x[1]].dc != True
        elif x[0] == "Line":
            return n.lines.loc[x[1]].dc != True
        else:
            logger.error("Unknown branch type for the instance {x}")

    return any(list(map(lambda x: is_ac_branch(x), ls)))


def simplify_links(n, costs, config, output, aggregation_strategies=dict()):
    ## Complex multi-node links are folded into end-points
    logger.info("Simplifying connected link components")

    if n.links.empty:
        with open(output.connection_costs, "w") as fp:
            pass
        return n, n.buses.index.to_series()

    # Determine connected link components, ignore all links but DC
    adjacency_matrix = n.adjacency_matrix(
        branch_components=["Link"],
        weights=dict(Link=(n.links.carrier == "DC").astype(float)),
    )

    _, labels = connected_components(adjacency_matrix, directed=False)
    labels = pd.Series(labels, n.buses.index)

    G = n.graph()

    # Split DC part by supernodes
    def split_links(nodes):
        nodes = frozenset(nodes)

        seen = set()
        supernodes = {m for m in nodes if len(G.adj[m]) > 2 or (set(G.adj[m]) - nodes)}

        for u in supernodes:
            for m, ls in G.adj[u].items():
                # AC lines can be captured in case of complicated network topologies
                # even despite using `nodes` defined by links
                if m not in nodes or m in seen or contains_ac(ls):
                    continue

                buses = [u, m]
                links = [list(ls)]  # [name for name in ls]]

                while m not in (supernodes | seen):
                    seen.add(m)
                    for m2, ls2 in G.adj[m].items():
                        # there may be AC lines which connect ends of DC chains
                        # TODO remove after debug
                        if m2 in seen or m2 == u or contains_ac(ls2):
                            continue
                        buses.append(m2)
                        links.append(list(ls2))  # [name for name in ls])
                        break
                    else:
                        # stub
                        break
                    m = m2
                if m != u:
                    yield pd.Index((u, m)), buses, links
            seen.add(u)

    busmap = n.buses.index.to_series()

    connection_costs_per_link = _prepare_connection_costs_per_link(n, costs, config)
    connection_costs_to_bus = pd.DataFrame(
        0.0, index=n.buses.index, columns=list(connection_costs_per_link)
    )

    for lbl in labels.value_counts().loc[lambda s: s > 2].index:
        for b, buses, links in split_links(labels.index[labels == lbl]):
            if len(buses) <= 2:
                continue

            logger.debug("nodes = {}".format(labels.index[labels == lbl]))
            logger.debug("b = {}\nbuses = {}\nlinks = {}".format(b, buses, links))

            m = sp.spatial.distance_matrix(
                n.buses.loc[b, ["x", "y"]], n.buses.loc[buses[1:-1], ["x", "y"]]
            )
            busmap.loc[buses] = b[np.r_[0, m.argmin(axis=0), 1]]
            connection_costs_to_bus.loc[buses] += _compute_connection_costs_to_bus(
                n, busmap, costs, config, connection_costs_per_link, buses
            )

            all_links = [i for _, i in sum(links, [])]

            p_max_pu = config["links"].get("p_max_pu", 1.0)
            lengths = n.links.loc[all_links, "length"]
            name = lengths.idxmax() + "+{}".format(len(links) - 1)
            params = dict(
                carrier="DC",
                bus0=b[0],
                bus1=b[1],
                length=sum(
                    n.links.loc[[i for _, i in l], "length"].mean() for l in links
                ),
                p_nom=min(n.links.loc[[i for _, i in l], "p_nom"].sum() for l in links),
                underwater_fraction=sum(
                    lengths
                    / lengths.sum()
                    * n.links.loc[all_links, "underwater_fraction"]
                ),
                p_max_pu=p_max_pu,
                p_min_pu=-p_max_pu,
                underground=False,
                under_construction=False,
            )

            logger.info(
                "Joining the links {} connecting the buses {} to simple link {}".format(
                    ", ".join(all_links), ", ".join(buses), name
                )
            )

            n.mremove("Link", all_links)

            static_attrs = n.components["Link"]["attrs"].loc[lambda df: df.static]
            for attr, default in static_attrs.default.items():
                params.setdefault(attr, default)
            n.links.loc[name] = params

    logger.debug("Collecting all components using the busmap")

    exclude_carriers = config["cluster_options"]["simplify_network"].get(
        "exclude_carriers", []
    )

    _aggregate_and_move_components(
        n,
        busmap,
        connection_costs_to_bus,
        output,
        aggregation_strategies=aggregation_strategies,
        exclude_carriers=exclude_carriers,
    )
    return n, busmap


def remove_stubs(n, costs, config, output, aggregation_strategies=dict()):
    logger.info("Removing stubs")

    across_borders = config["cluster_options"]["simplify_network"].get(
        "remove_stubs_across_borders", True
    )
    matching_attrs = [] if across_borders else ["country"]

    busmap = busmap_by_stubs(n, matching_attrs)

    connection_costs_to_bus = _compute_connection_costs_to_bus(n, busmap, costs, config)

    exclude_carriers = config["cluster_options"]["simplify_network"].get(
        "exclude_carriers", []
    )

    _aggregate_and_move_components(
        n,
        busmap,
        connection_costs_to_bus,
        output,
        aggregation_strategies=aggregation_strategies,
        exclude_carriers=exclude_carriers,
    )

    return n, busmap


def aggregate_to_substations(n, aggregation_strategies=dict(), buses_i=None):
    # can be used to aggregate a selection of buses to electrically closest neighbors
    # if no buses are given, nodes that are no substations or without offshore connection are aggregated

    if buses_i is None:
        logger.info(
            "Aggregating buses that are no substations or have no valid offshore connection"
        )
        # isolated buses should be excluded from adjacency_matrix calculations
        n.determine_network_topology()
        # duplicated sub-networks mean that there is at least one interconnection between buses
        i_islands = n.buses[
            (~n.buses.duplicated(subset=["sub_network"], keep=False))
            & (n.buses.carrier == "AC")
        ].index

        buses_i = list(
            set(n.buses.index)
            - set(i_islands)
            - set(n.generators.bus)
            - set(n.loads.bus)
        )

    weight = pd.concat(
        {
            "Line": n.lines.length / n.lines.s_nom.clip(1e-3),
            "Link": n.links.length / n.links.p_nom.clip(1e-3),
        }
    )

    adj = n.adjacency_matrix(branch_components=["Line", "Link"], weights=weight)

    bus_indexer = n.buses.index.get_indexer(buses_i)
    dist = pd.DataFrame(
        dijkstra(adj, directed=False, indices=bus_indexer), buses_i, n.buses.index
    )

    dist[
        buses_i
    ] = np.inf  # bus in buses_i should not be assigned to different bus in buses_i

    # avoid assignment a bus to a wrong country
    for c in n.buses.country.unique():
        incountry_b = n.buses.country == c
        dist.loc[incountry_b, ~incountry_b] = np.inf

    # avoid assignment DC buses to AC ones
    for c in n.buses.carrier.unique():
        incarrier_b = n.buses.carrier == c
        dist.loc[incarrier_b, ~incarrier_b] = np.inf

    busmap = n.buses.index.to_series()
    if not dist.empty:
        busmap.loc[buses_i] = dist.idxmin(1)

    bus_strategies, generator_strategies = get_aggregation_strategies(
        aggregation_strategies
    )

    clustering = get_clustering_from_busmap(
        n,
        busmap,
        bus_strategies=bus_strategies,
        aggregate_generators_weighted=True,
        aggregate_generators_carriers=None,
        aggregate_one_ports=["Load", "StorageUnit"],
        line_length_factor=1.0,
        generator_strategies=generator_strategies,
        scale_link_capital_costs=False,
    )
    return clustering.network, busmap


def cluster(
    n, n_clusters, config, algorithm="hac", feature=None, aggregation_strategies=dict()
):
    logger.info(f"Clustering to {n_clusters} buses")

    focus_weights = config.get("focus_weights", None)
    alternative_clustering = config["cluster_options"]["alternative_clustering"]
    gadm_layer_id = config["build_shape_options"]["gadm_layer_id"]
    geo_crs = config["crs"]["geo_crs"]
    country_list = config["countries"]

    renewable_carriers = pd.Index(
        [
            tech
            for tech in n.generators.carrier.unique()
            if tech.split("-", 2)[0] in config["renewable"]
        ]
    )

    def consense(x):
        v = x.iat[0]
        assert (
            x == v
        ).all() or x.isnull().all(), "The `potential` configuration option must agree for all renewable carriers, for now!"
        return v

    potential_mode = (
        consense(
            pd.Series(
                [config["renewable"][tech]["potential"] for tech in renewable_carriers]
            )
        )
        if len(renewable_carriers) > 0
        else "conservative"
    )
    clustering = clustering_for_n_clusters(
        n,
        n_clusters,
        alternative_clustering,
        gadm_layer_id,
        geo_crs,
        country_list,
        custom_busmap=False,
        aggregation_strategies=aggregation_strategies,
        potential_mode=potential_mode,
        solver_name=config["solving"]["solver"]["name"],
        algorithm=algorithm,
        feature=feature,
        focus_weights=focus_weights,
    )

    return clustering.network, clustering.busmap


def drop_isolated_nodes(n, threshold):
    """
    Find isolated nodes in the network and drop those of them which don't have
    load or have a load value lower than a threshold.

    Parameters
    ----------
    iso_code : PyPSA.Network
        Original network
    threshold: float
        Load power used as a threshold to drop isolated nodes

    Returns
    -------
    modified network
    """
    n.determine_network_topology()

    # keep original values of the overall load and generation in the network
    # to track changes due to drop of buses
    generators_mean_origin = n.generators.p_nom.mean()
    load_mean_origin = n.loads_t.p_set.mean().mean()

    # duplicated sub-networks mean that there is at least one interconnection between buses
    off_buses = n.buses[
        (~n.buses.duplicated(subset=["sub_network"], keep=False))
        & (n.buses.carrier == "AC")
    ]
    i_islands = off_buses.index

    if off_buses.empty:
        return n

    # skip a isolated bus for countries represented by only isolated nodes
    for c in off_buses["country"].unique():
        isc_offbus = off_buses["country"] == c
        n_bus_off = isc_offbus.sum()
        n_bus = (n.buses["country"] == c).sum()
        if n_bus_off == n_bus:
            i_islands = i_islands[i_islands != off_buses[isc_offbus].index[0]]

    i_load_islands = n.loads_t.p_set.columns.intersection(i_islands)

    # isolated buses without load should be discarded
    isl_no_load = i_islands.difference(i_load_islands)

    # isolated buses with load lower than a specified threshold should be discarded as well
    i_small_load = i_load_islands[
        n.loads_t.p_set[i_load_islands].mean(axis=0) <= threshold
    ]

    if (i_islands.empty) | (len(i_small_load) == 0):
        return n

    i_to_drop = isl_no_load.to_list() + i_small_load.to_list()

    i_loads_drop = n.loads[n.loads.bus.isin(i_to_drop)].index
    i_generators_drop = n.generators[n.generators.bus.isin(i_to_drop)].index

    n.mremove("Bus", i_to_drop)
    n.mremove("Load", i_loads_drop)
    n.mremove("Generator", i_generators_drop)

    n.determine_network_topology()

    load_mean_final = n.loads_t.p_set.mean().mean()
    generators_mean_final = n.generators.p_nom.mean()

    logger.info(
        f"Dropped {len(i_to_drop)} buses. A resulted load discrepancy is {(100 * ((load_mean_final - load_mean_origin)/load_mean_origin)):2.1}% and {(100 * ((generators_mean_final - generators_mean_origin)/generators_mean_origin)):2.1}% for average load and generation capacity, respectively"
    )

    return n


def merge_isolated_nodes(n, threshold, aggregation_strategies=dict()):
    """
    Find isolated nodes in the network and merge those of them which have load
    value below than a specified threshold into a single isolated node which
    represents all the remote generation.

    Parameters
    ----------
    n : PyPSA.Network
        Original network
    threshold : float
        Load power used as a threshold to merge isolated nodes

    Returns
    -------
    modified network
    """
    # keep original values of the overall load and generation in the network
    # to track changes due to drop of buses
    generators_mean_origin = n.generators.p_nom.mean()
    load_mean_origin = n.loads_t.p_set.mean().mean()

    n.determine_network_topology()

    # duplicated sub-networks mean that there is at least one interconnection between buses
    i_islands = n.buses[
        (~n.buses.duplicated(subset=["sub_network"], keep=False))
        & (n.buses.carrier == "AC")
    ].index

    # isolated buses with load below than a specified threshold should be merged
    i_load_islands = n.loads_t.p_set.columns.intersection(i_islands)
    i_suffic_load = i_load_islands[
        n.loads_t.p_set[i_load_islands].mean(axis=0) <= threshold
    ]

    # all the nodes to be merged should be mapped into a single node
    map_isolated_node_by_country = (
        n.buses.assign(bus_id=n.buses.index)
        .loc[i_suffic_load]
        .groupby("country")["bus_id"]
        .first()
        .to_dict()
    )
    isolated_buses_mapping = n.buses.loc[i_suffic_load, "country"].replace(
        map_isolated_node_by_country
    )
    busmap = (
        n.buses.index.to_series()
        .replace(isolated_buses_mapping)
        .astype(str)
        .rename("busmap")
    )

    # return the original network if no changes are detected
    if (busmap.index == busmap).all():
        return n, n.buses.index.to_series()

    bus_strategies, generator_strategies = get_aggregation_strategies(
        aggregation_strategies
    )

    clustering = get_clustering_from_busmap(
        n,
        busmap,
        bus_strategies=bus_strategies,
        aggregate_generators_weighted=True,
        aggregate_generators_carriers=None,
        aggregate_one_ports=["Load", "StorageUnit"],
        line_length_factor=1.0,
        generator_strategies=generator_strategies,
        scale_link_capital_costs=False,
    )

    load_mean_final = n.loads_t.p_set.mean().mean()
    generators_mean_final = n.generators.p_nom.mean()

    logger.info(
        f"Merged {len(i_suffic_load)} buses. Load attached to a single bus with discrepancies of {(100 * ((load_mean_final - load_mean_origin)/load_mean_origin)):2.1E}% and {(100 * ((generators_mean_final - generators_mean_origin)/generators_mean_origin)):2.1E}% for load and generation capacity, respectively"
    )

    return clustering.network, busmap


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("simplify_network", simpl="")
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)

    base_voltage = snakemake.config["electricity"]["base_voltage"]
    linetype = snakemake.config["lines"]["types"][base_voltage]

    aggregation_strategies = snakemake.config["cluster_options"].get(
        "aggregation_strategies", {}
    )
    # translate str entries of aggregation_strategies to pd.Series functions:
    aggregation_strategies = {
        p: {k: getattr(pd.Series, v) for k, v in aggregation_strategies[p].items()}
        for p in aggregation_strategies.keys()
    }
    n, trafo_map = simplify_network_to_base_voltage(n, linetype, base_voltage)

    Nyears = n.snapshot_weightings.objective.sum() / 8760

    technology_costs = load_costs(
        snakemake.input.tech_costs,
        snakemake.config["costs"],
        snakemake.config["electricity"],
        Nyears,
    )

    n, simplify_links_map = simplify_links(
        n, technology_costs, snakemake.config, snakemake.output, aggregation_strategies
    )

    busmaps = [trafo_map, simplify_links_map]

    cluster_config = snakemake.config["cluster_options"]["simplify_network"]
    if cluster_config.get("remove_stubs", True):
        n, stub_map = remove_stubs(
            n,
            technology_costs,
            snakemake.config,
            snakemake.output,
            aggregation_strategies=aggregation_strategies,
        )
        busmaps.append(stub_map)

    if cluster_config.get("to_substations", False):
        n, substation_map = aggregate_to_substations(n, aggregation_strategies)
        busmaps.append(substation_map)

    # treatment of outliers (nodes without a profile for considered carrier):
    # all nodes that have no profile of the given carrier are being aggregated to closest neighbor
    if (
        snakemake.config.get("cluster_options", {})
        .get("cluster_network", {})
        .get("algorithm", "hac")
        == "hac"
        or cluster_config.get("algorithm", "hac") == "hac"
    ):
        carriers = (
            cluster_config.get("feature", "solar+onwind-time").split("-")[0].split("+")
        )
        for carrier in carriers:
            # isolated buses should be excluded from adjacency_matrix calculations
            n.determine_network_topology()
            # duplicated sub-networks mean that there is at least one interconnection between buses
            i_islands = n.buses[
                (~n.buses.duplicated(subset=["sub_network"], keep=False))
                & (n.buses.carrier == "AC")
            ].index

            buses_i = list(
                set(n.buses.index)
                - set(i_islands)
                - set(n.generators.query("carrier == @carrier").bus)
            )
            logger.info(
                f"clustering preparation (hac): aggregating {len(buses_i)} buses of type {carrier}."
            )
            n, busmap_hac = aggregate_to_substations(n, aggregation_strategies, buses_i)
            busmaps.append(busmap_hac)

    if snakemake.wildcards.simpl:
        n, cluster_map = cluster(
            n,
            int(snakemake.wildcards.simpl),
            snakemake.config,
            cluster_config.get("algorithm", "hac"),
            cluster_config.get("feature", None),
            aggregation_strategies,
        )
        busmaps.append(cluster_map)

    # some entries in n.buses are not updated in previous functions, therefore can be wrong. as they are not needed
    # and are lost when clustering (for example with the simpl wildcard), we remove them for consistency:
    buses_c = {
        "symbol",
        "tags",
        "under_construction",
        "substation_lv",
        "substation_off",
    }.intersection(n.buses.columns)
    n.buses = n.buses.drop(buses_c, axis=1)

    update_p_nom_max(n)

    p_threshold_drop_isolated = max(
        0.0, cluster_config.get("p_threshold_drop_isolated", 0.0)
    )
    p_threshold_merge_isolated = cluster_config.get("p_threshold_merge_isolated", False)

    n = drop_isolated_nodes(n, threshold=p_threshold_drop_isolated)
    if p_threshold_merge_isolated:
        n, merged_nodes_map = merge_isolated_nodes(
            n,
            threshold=p_threshold_merge_isolated,
            aggregation_strategies=aggregation_strategies,
        )
        busmaps.append(merged_nodes_map)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output.network)

    busmap_s = reduce(lambda x, y: x.map(y), busmaps[1:], busmaps[0])
    busmap_s.to_csv(snakemake.output.busmap)

    cluster_regions(busmaps, snakemake.input, snakemake.output)
