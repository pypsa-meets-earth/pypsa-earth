# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors, 2021 PyPSA-Africa Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later
# coding: utf-8
"""
Creates networks clustered to ``{cluster}`` number of zones with aggregated buses, generators and transmission corridors.

Relevant Settings
-----------------

.. code:: yaml

    clustering:
        aggregation_strategies:

    focus_weights:

    solving:
        solver:
            name:

    lines:
        length_factor:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`toplevel_cf`, :ref:`renewable_cf`, :ref:`solving_cf`, :ref:`lines_cf`

Inputs
------

- ``resources/regions_onshore_elec_s{simpl}.geojson``: confer :ref:`simplify`
- ``resources/regions_offshore_elec_s{simpl}.geojson``: confer :ref:`simplify`
- ``resources/busmap_elec_s{simpl}.csv``: confer :ref:`simplify`
- ``networks/elec_s{simpl}.nc``: confer :ref:`simplify`
- ``data/custom_busmap_elec_s{simpl}_{clusters}.csv``: optional input

Outputs
-------

- ``resources/regions_onshore_elec_s{simpl}_{clusters}.geojson``:

    .. image:: ../img/regions_onshore_elec_s_X.png
        :scale: 33 %

- ``resources/regions_offshore_elec_s{simpl}_{clusters}.geojson``:

    .. image:: ../img/regions_offshore_elec_s_X.png
        :scale: 33 %

- ``resources/busmap_elec_s{simpl}_{clusters}.csv``: Mapping of buses from ``networks/elec_s{simpl}.nc`` to ``networks/elec_s{simpl}_{clusters}.nc``;
- ``resources/linemap_elec_s{simpl}_{clusters}.csv``: Mapping of lines from ``networks/elec_s{simpl}.nc`` to ``networks/elec_s{simpl}_{clusters}.nc``;
- ``networks/elec_s{simpl}_{clusters}.nc``:

    .. image:: ../img/elec_s_X.png
        :scale: 40  %

Description
-----------

.. note::

    **Why is clustering used both in** ``simplify_network`` **and** ``cluster_network`` **?**

        Consider for example a network ``networks/elec_s100_50.nc`` in which
        ``simplify_network`` clusters the network to 100 buses and in a second
        step ``cluster_network``` reduces it down to 50 buses.

        In preliminary tests, it turns out, that the principal effect of
        changing spatial resolution is actually only partially due to the
        transmission network. It is more important to differentiate between
        wind generators with higher capacity factors from those with lower
        capacity factors, i.e. to have a higher spatial resolution in the
        renewable generation than in the number of buses.

        The two-step clustering allows to study this effect by looking at
        networks like ``networks/elec_s100_50m.nc``. Note the additional
        ``m`` in the ``{cluster}`` wildcard. So in the example network
        there are still up to 100 different wind generators.

        In combination these two features allow you to study the spatial
        resolution of the transmission network separately from the
        spatial resolution of renewable generators.

    **Is it possible to run the model without the** ``simplify_network`` **rule?**

        No, the network clustering methods in the PyPSA module
        `pypsa.networkclustering <https://github.com/PyPSA/PyPSA/blob/master/pypsa/networkclustering.py>`_
        do not work reliably with multiple voltage levels and transformers.

.. tip::
    The rule :mod:`cluster_all_networks` runs
    for all ``scenario`` s in the configuration file
    the rule :mod:`cluster_network`.

Exemplary unsolved network clustered to 512 nodes:

.. image:: ../img/elec_s_512.png
    :scale: 40  %
    :align: center

Exemplary unsolved network clustered to 256 nodes:

.. image:: ../img/elec_s_256.png
    :scale: 40  %
    :align: center

Exemplary unsolved network clustered to 128 nodes:

.. image:: ../img/elec_s_128.png
    :scale: 40  %
    :align: center

Exemplary unsolved network clustered to 37 nodes:

.. image:: ../img/elec_s_37.png
    :scale: 40  %
    :align: center

"""
import logging
import os
from functools import reduce

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyomo.environ as po
import pypsa
import seaborn as sns
import shapely
from _helpers import (
    REGION_COLS,
    configure_logging,
    get_aggregation_strategies,
    sets_path_to_root,
    two_2_three_digits_country,
    update_p_nom_max,
)
from add_electricity import load_costs
from build_shapes import add_gdp_data, add_population_data, get_GADM_layer
from pypsa.networkclustering import (
    _make_consense,
    busmap_by_kmeans,
    get_clustering_from_busmap,
)
from shapely.geometry import Point

idx = pd.IndexSlice

logger = logging.getLogger(__name__)


def normed(x):
    return (x / x.sum()).fillna(0.0)


def weighting_for_country(n, x):
    conv_carriers = {"OCGT", "CCGT", "PHS", "hydro"}
    gen = n.generators.loc[n.generators.carrier.isin(conv_carriers)].groupby(
        "bus"
    ).p_nom.sum().reindex(n.buses.index, fill_value=0.0) + n.storage_units.loc[
        n.storage_units.carrier.isin(conv_carriers)
    ].groupby(
        "bus"
    ).p_nom.sum().reindex(
        n.buses.index, fill_value=0.0
    )
    load = n.loads_t.p_set.mean().groupby(n.loads.bus).sum()

    b_i = x.index
    g = normed(gen.reindex(b_i, fill_value=0))
    l = normed(load.reindex(b_i, fill_value=0))

    w = g + l

    if w.max() == 0.0:
        logger.warning(
            f"Null weighting for buses of country {x.country.iloc[0]}: returned default uniform weighting"
        )
        return pd.Series(1.0, index=w.index)
    else:
        return (w * (100.0 / w.max())).clip(lower=1.0).astype(int)


def distribute_clusters(
    inputs, config, n, n_clusters, focus_weights=None, solver_name=None
):
    """Determine the number of clusters per country"""

    distribution_cluster = config["cluster_options"]["distribute_cluster"]
    country_list = config["countries"]
    year = config["build_shape_options"]["year"]
    update = config["build_shape_options"]["update_file"]
    out_logging = config["build_shape_options"]["out_logging"]
    nprocesses = config["build_shape_options"]["nprocesses"]

    if solver_name is None:
        solver_name = config["solving"]["solver"]["name"]

    if distribution_cluster == ["load"]:
        L = (
            n.loads_t.p_set.mean()
            .groupby(n.loads.bus)
            .sum()
            .groupby([n.buses.country])
            .sum()
            .pipe(normed)
        )
        assert len(L.index) == len(n.buses.country.unique()), (
            "The following countries have no load: "
            f"{list(set(L.index).symmetric_difference(set(n.buses.country.unique())))}"
        )
        distribution_factor = L

    if distribution_cluster == ["pop"]:
        df_pop_c = gpd.read_file(inputs.country_shapes).rename(
            columns={"name": "country"}
        )
        add_population_data(
            df_pop_c, country_list, year, update, out_logging, nprocesses=nprocesses
        )
        P = df_pop_c.loc[:, ("country", "pop")]
        P = P.groupby(P["country"]).sum().pipe(normed).squeeze()
        distribution_factor = P

    if distribution_cluster == ["gdp"]:
        df_gdp_c = gpd.read_file(inputs.country_shapes).rename(
            columns={"name": "country"}
        )
        add_gdp_data(
            df_gdp_c,
            year,
            update,
            out_logging,
            name_file_nc="GDP_PPP_1990_2015_5arcmin_v2.nc",
        )
        G = df_gdp_c.loc[:, ("country", "gdp")]
        G = G.groupby(df_gdp_c["country"]).sum().pipe(normed).squeeze()
        distribution_factor = G

    # TODO: 1. Check if sub_networks can be added here i.e. ["country", "sub_network"]
    N = n.buses.groupby(["country"]).size()

    assert (
        n_clusters >= len(N) and n_clusters <= N.sum()
    ), f"Number of clusters must be {len(N)} <= n_clusters <= {N.sum()} for this selection of countries."

    if focus_weights is not None:
        total_focus = sum(list(focus_weights.values()))

        assert (
            total_focus <= 1.0
        ), "The sum of focus weights must be less than or equal to 1."

        for country, weight in focus_weights.items():
            distribution_factor[country] = weight / len(distribution_factor[country])

        remainder = [
            c not in focus_weights.keys()
            for c in distribution_factor.index.get_level_values("country")
        ]
        distribution_factor[remainder] = distribution_factor.loc[remainder].pipe(
            normed
        ) * (1 - total_focus)

        logger.warning("Using custom focus weights for determining number of clusters.")

    assert np.isclose(
        distribution_factor.sum(), 1.0, rtol=1e-3
    ), f"Country weights L must sum up to 1.0 when distributing clusters. Is {distribution_factor.sum()}."

    m = po.ConcreteModel()

    def n_bounds(model, n_id):
        """
        Create a function that makes a bound pair for pyomo

        Use n_bounds(model, n_id) if N is Single-Index
        Use n_bounds(model, *n_id) if N is Multi-Index
        Example: https://pyomo.readthedocs.io/en/stable/pyomo_modeling_components/Variables.html

        Returns
        -------
        bounds = A function (or Python object) that gives a (lower,upper) bound pair i.e.(1,10) for the variable
        """
        return (1, N[n_id])

    m.n = po.Var(list(distribution_factor.index), bounds=n_bounds, domain=po.Integers)
    m.tot = po.Constraint(expr=(po.summation(m.n) == n_clusters))
    m.objective = po.Objective(
        expr=sum(
            (m.n[i] - distribution_factor.loc[i] * n_clusters) ** 2
            for i in distribution_factor.index
        ),
        sense=po.minimize,
    )

    opt = po.SolverFactory(solver_name)
    if not opt.has_capability("quadratic_objective"):
        logger.warning(
            f"The configured solver `{solver_name}` does not support quadratic objectives. Falling back to `ipopt`."
        )
        opt = po.SolverFactory("ipopt")

    results = opt.solve(m)
    assert (
        results["Solver"][0]["Status"] == "ok"
    ), f"Solver returned non-optimally: {results}"

    return (
        pd.Series(m.n.get_values(), index=distribution_factor.index).round().astype(int)
    )


def busmap_for_gadm_clusters(inputs, n, gadm_level, geo_crs, country_list):
    # gdf = get_GADM_layer(country_list, gadm_level, geo_crs)
    gdf = gpd.read_file(inputs.gadm_shapes)

    def locate_bus(coords, co):
        gdf_co = gdf[gdf["GADM_ID"].str.contains(two_2_three_digits_country(co))]
        point = Point(coords["x"], coords["y"])

        try:
            return gdf_co[gdf_co.contains(point)]["GADM_ID"].item()

        except ValueError:
            return gdf_co[
                gdf_co.geometry == min(gdf_co.geometry, key=(point.distance))
            ]["GADM_ID"].item()

    buses = n.buses
    buses["gadm_{}".format(gadm_level)] = buses[["x", "y", "country"]].apply(
        lambda bus: locate_bus(bus[["x", "y"]], bus["country"]), axis=1
    )
    busmap = buses["gadm_{}".format(gadm_level)]

    return busmap


def busmap_for_n_clusters(
    inputs,
    config,
    n,
    n_clusters,
    solver_name,
    focus_weights=None,
    algorithm="kmeans",
    **algorithm_kwds,
):
    if algorithm == "kmeans":
        algorithm_kwds.setdefault("n_init", 1000)
        algorithm_kwds.setdefault("max_iter", 30000)
        algorithm_kwds.setdefault("tol", 1e-6)
        algorithm_kwds.setdefault("random_state", 0)

    n.determine_network_topology()
    # n.lines.loc[:, "sub_network"] = "0"  # current fix

    if n.buses.country.nunique() > 1:
        n_clusters = distribute_clusters(
            inputs,
            config,
            n,
            n_clusters,
            focus_weights=focus_weights,
            solver_name=solver_name,
        )

    def reduce_network(n, buses):
        nr = pypsa.Network()
        nr.import_components_from_dataframe(buses, "Bus")
        nr.import_components_from_dataframe(
            n.lines.loc[
                n.lines.bus0.isin(buses.index) & n.lines.bus1.isin(buses.index)
            ],
            "Line",
        )
        return nr

    def busmap_for_country(x):
        # A number of the countries in the clustering can be > 1
        if isinstance(n_clusters, pd.Series):
            if isinstance(x.name, tuple):
                n_cluster_c = n_clusters[x.name[0]]
                prefix = x.name[0] + x.name[1] + " "
            else:
                n_cluster_c = n_clusters[x.name]
                prefix = x.name + " "
        else:
            n_cluster_c = n_clusters
            prefix = x.name[0] + x.name[1] + " "

        logger.debug(f"Determining busmap for country {prefix[:-1]}")
        if len(x) == 1:
            return pd.Series(prefix + "0", index=x.index)
        weight = weighting_for_country(n, x)

        if algorithm == "kmeans":
            return prefix + busmap_by_kmeans(
                n, weight, n_cluster_c, buses_i=x.index, **algorithm_kwds
            )

        else:
            raise ValueError(
                f"`algorithm` must be one of 'kmeans' or 'louvain'. Is {algorithm}."
            )

    return (
        n.buses.groupby(
            ["country"],
            # ["country", "sub_network"] # TODO: 2. Add sub_networks (see previous TODO)
            group_keys=False,
        )
        .apply(busmap_for_country)
        .squeeze(axis=0)
        .rename("busmap")
    )


def clustering_for_n_clusters(
    inputs,
    config,
    n,
    n_clusters,
    alternative_clustering,
    gadm_layer_id,
    geo_crs,
    country_list,
    custom_busmap=False,
    aggregate_carriers=None,
    line_length_factor=1.25,
    aggregation_strategies=dict(),
    solver_name="cbc",
    algorithm="kmeans",
    extended_link_costs=0,
    focus_weights=None,
):
    bus_strategies, generator_strategies = get_aggregation_strategies(
        aggregation_strategies
    )

    if not isinstance(custom_busmap, pd.Series):
        if alternative_clustering:
            busmap = busmap_for_gadm_clusters(
                inputs, n, gadm_layer_id, geo_crs, country_list
            )
        else:
            busmap = busmap_for_n_clusters(
                inputs, config, n, n_clusters, solver_name, focus_weights, algorithm
            )
    else:
        busmap = custom_busmap

    clustering = get_clustering_from_busmap(
        n,
        busmap,
        bus_strategies=bus_strategies,
        aggregate_generators_weighted=True,
        aggregate_generators_carriers=aggregate_carriers,
        aggregate_one_ports=["Load", "StorageUnit"],
        line_length_factor=line_length_factor,
        generator_strategies=generator_strategies,
        scale_link_capital_costs=False,
    )

    if not n.links.empty:
        nc = clustering.network
        nc.links["underwater_fraction"] = (
            n.links.eval("underwater_fraction * length").div(nc.links.length).dropna()
        )
        nc.links["capital_cost"] = nc.links["capital_cost"].add(
            (nc.links.length - n.links.length).clip(lower=0).mul(extended_link_costs),
            fill_value=0,
        )

    return clustering


def save_to_geojson(s, fn):
    if os.path.exists(fn):
        os.unlink(fn)
    df = s.reset_index()
    schema = {**gpd.io.file.infer_schema(df), "geometry": "Unknown"}
    df.to_file(fn, driver="GeoJSON", schema=schema)


def cluster_regions(busmaps, inputs, output):
    busmap = reduce(lambda x, y: x.map(y), busmaps[1:], busmaps[0])

    for which in ("regions_onshore", "regions_offshore"):
        # regions = gpd.read_file(getattr(input, which)).set_index("name")
        regions = gpd.read_file(getattr(inputs, which))
        regions = regions.reindex(columns=REGION_COLS).set_index("name")
        aggfunc = dict(x="mean", y="mean", country="first")
        regions_c = regions.dissolve(busmap, aggfunc=aggfunc)
        regions_c.index.name = "name"
        regions_c = regions_c.reset_index()
        regions_c.to_file(getattr(output, which))


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "cluster_network", network="elec", simpl="", clusters="60"
        )
        sets_path_to_root("pypsa-earth")
    configure_logging(snakemake)

    inputs, outputs, config = snakemake.input, snakemake.output, snakemake.config

    n = pypsa.Network(inputs.network)

    alternative_clustering = snakemake.config["cluster_options"][
        "alternative_clustering"
    ]
    gadm_layer_id = snakemake.config["build_shape_options"]["gadm_layer_id"]
    focus_weights = snakemake.config.get("focus_weights", None)
    country_list = snakemake.config["countries"]
    geo_crs = snakemake.config["crs"]["geo_crs"]

    renewable_carriers = pd.Index(
        [
            tech
            for tech in n.generators.carrier.unique()
            if tech in snakemake.config["renewable"]  # TODO ror not cap
        ]
    )

    if snakemake.wildcards.clusters.endswith("m"):
        n_clusters = int(snakemake.wildcards.clusters[:-1])
        aggregate_carriers = snakemake.config["electricity"].get(
            "conventional_carriers"
        )
    else:
        n_clusters = int(snakemake.wildcards.clusters)
        aggregate_carriers = None

    if n_clusters == len(n.buses):
        # Fast-path if no clustering is necessary
        busmap = n.buses.index.to_series()
        linemap = n.lines.index.to_series()
        clustering = pypsa.networkclustering.Clustering(
            n, busmap, linemap, linemap, pd.Series(dtype="O")
        )
    elif len(n.buses) < n_clusters:
        logger.error(
            f"Desired number of clusters ({n_clusters}) higher than the number of buses ({len(n.buses)})"
        )
    else:
        line_length_factor = snakemake.config["lines"]["length_factor"]
        Nyears = n.snapshot_weightings.objective.sum() / 8760
        hvac_overhead_cost = load_costs(
            snakemake.input.tech_costs,
            snakemake.config["costs"],
            snakemake.config["electricity"],
            Nyears,
        ).at["HVAC overhead", "capital_cost"]

        def consense(x):
            v = x.iat[0]
            assert (
                x == v
            ).all() or x.isnull().all(), "The `potential` configuration option must agree for all renewable carriers, for now!"
            return v

        aggregation_strategies = snakemake.config["cluster_options"].get(
            "aggregation_strategies", {}
        )
        # translate str entries of aggregation_strategies to pd.Series functions:
        aggregation_strategies = {
            p: {k: getattr(pd.Series, v) for k, v in aggregation_strategies[p].items()}
            for p in aggregation_strategies.keys()
        }
        custom_busmap = snakemake.config["enable"].get("custom_busmap", False)
        if custom_busmap:
            busmap = pd.read_csv(
                snakemake.input.custom_busmap, index_col=0, squeeze=True
            )
            busmap.index = busmap.index.astype(str)
            logger.info(f"Imported custom busmap from {snakemake.input.custom_busmap}")
        cluster_config = snakemake.config.get("cluster_options", {}).get(
            "cluster_network", {}
        )
        clustering = clustering_for_n_clusters(
            inputs,
            config,
            n,
            n_clusters,
            alternative_clustering,
            gadm_layer_id,
            geo_crs,
            country_list,
            custom_busmap,
            aggregate_carriers,
            line_length_factor,
            aggregation_strategies,
            solver_name=snakemake.config["solving"]["solver"]["name"],
            extended_link_costs=hvac_overhead_cost,
            focus_weights=focus_weights,
        )

    update_p_nom_max(clustering.network)
    clustering.network.meta = dict(
        snakemake.config, **dict(wildcards=dict(snakemake.wildcards))
    )
    clustering.network.export_to_netcdf(outputs.network)
    for attr in (
        "busmap",
        "linemap",
    ):  # also available: linemap_positive, linemap_negative
        getattr(clustering, attr).to_csv(outputs[attr])

    cluster_regions((clustering.busmap,), inputs, outputs)
