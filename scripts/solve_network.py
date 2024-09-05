# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Solves linear optimal power flow for a network iteratively while updating
reactances.

Relevant Settings
-----------------

.. code:: yaml

    solving:
        tmpdir:
        options:
            formulation:
            clip_p_max_pu:
            load_shedding:
            noisy_costs:
            nhours:
            min_iterations:
            max_iterations:
            skip_iterations:
            track_iterations:
        solver:
            name:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`electricity_cf`, :ref:`solving_cf`, :ref:`plotting_cf`

Inputs
------

- ``networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc``: confer :ref:`prepare`

Outputs
-------

- ``results/networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc``: Solved PyPSA network including optimisation results

    .. image:: /img/results.png
        :width: 40 %

Description
-----------

Total annual system costs are minimised with PyPSA. The full formulation of the
linear optimal power flow (plus investment planning)
is provided in the
`documentation of PyPSA <https://pypsa.readthedocs.io/en/latest/optimal_power_flow.html#linear-optimal-power-flow>`_.
The optimization is based on the ``pyomo=False`` setting in the :func:`network.lopf` and  :func:`pypsa.linopf.ilopf` function.
Additionally, some extra constraints specified in :mod:`prepare_network` are added.

Solving the network in multiple iterations is motivated through the dependence of transmission line capacities and impedances on values of corresponding flows.
As lines are expanded their electrical parameters change, which renders the optimisation bilinear even if the power flow
equations are linearized.
To retain the computational advantage of continuous linear programming, a sequential linear programming technique
is used, where in between iterations the line impedances are updated.
Details (and errors made through this heuristic) are discussed in the paper

- Fabian Neumann and Tom Brown. `Heuristics for Transmission Expansion Planning in Low-Carbon Energy System Models <https://arxiv.org/abs/1907.10548>`_), *16th International Conference on the European Energy Market*, 2019. `arXiv:1907.10548 <https://arxiv.org/abs/1907.10548>`_.

.. warning::
    Capital costs of existing network components are not included in the objective function,
    since for the optimisation problem they are just a constant term (no influence on optimal result).

    Therefore, these capital costs are not included in ``network.objective``!

    If you want to calculate the full total annual system costs add these to the objective value.

.. tip::
    The rule :mod:`solve_all_networks` runs
    for all ``scenario`` s in the configuration file
    the rule :mod:`solve_network`.
"""
import logging
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
from _helpers import configure_logging, create_logger, override_component_attrs
from pypsa.descriptors import get_switchable_as_dense as get_as_dense
from pypsa.linopf import (
    define_constraints,
    define_variables,
    get_var,
    ilopf,
    join_exprs,
    linexpr,
    network_lopf,
)
from pypsa.linopt import define_constraints, get_var, join_exprs, linexpr

logger = create_logger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)


def prepare_network(n, solve_opts):
    if "clip_p_max_pu" in solve_opts:
        for df in (
            n.generators_t.p_max_pu,
            n.generators_t.p_min_pu,
            n.storage_units_t.inflow,
        ):
            df.where(df > solve_opts["clip_p_max_pu"], other=0.0, inplace=True)

    if "lv_limit" in n.global_constraints.index:
        n.line_volume_limit = n.global_constraints.at["lv_limit", "constant"]
        n.line_volume_limit_dual = n.global_constraints.at["lv_limit", "mu"]

    if solve_opts.get("load_shedding"):
        n.add("Carrier", "Load")
        n.madd(
            "Generator",
            n.buses.index,
            " load",
            bus=n.buses.index,
            carrier="load",
            sign=1e-3,  # Adjust sign to measure p and p_nom in kW instead of MW
            marginal_cost=1e2,  # Eur/kWh
            # intersect between macroeconomic and surveybased
            # willingness to pay
            # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
            p_nom=1e9,  # kW
        )

    if solve_opts.get("noisy_costs"):
        for t in n.iterate_components():
            # if 'capital_cost' in t.df:
            #    t.df['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.df)) - 0.5)
            if "marginal_cost" in t.df:
                np.random.seed(174)
                t.df["marginal_cost"] += 1e-2 + 2e-3 * (
                    np.random.random(len(t.df)) - 0.5
                )

        for t in n.iterate_components(["Line", "Link"]):
            np.random.seed(123)
            t.df["capital_cost"] += (
                1e-1 + 2e-2 * (np.random.random(len(t.df)) - 0.5)
            ) * t.df["length"]

    if solve_opts.get("nhours"):
        nhours = solve_opts["nhours"]
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760.0 / nhours

    if snakemake.config["foresight"] == "myopic":
        add_land_use_constraint(n)

    return n


def add_CCL_constraints(n, config):
    agg_p_nom_limits = config["electricity"].get("agg_p_nom_limits")

    try:
        agg_p_nom_minmax = pd.read_csv(agg_p_nom_limits, index_col=list(range(2)))
    except IOError:
        logger.exception(
            "Need to specify the path to a .csv file containing "
            "aggregate capacity limits per country in "
            "config['electricity']['agg_p_nom_limit']."
        )
    logger.info(
        "Adding per carrier generation capacity constraints for " "individual countries"
    )

    gen_country = n.generators.bus.map(n.buses.country)
    # cc means country and carrier
    p_nom_per_cc = (
        pd.DataFrame(
            {
                "p_nom": linexpr((1, get_var(n, "Generator", "p_nom"))),
                "country": gen_country,
                "carrier": n.generators.carrier,
            }
        )
        .dropna(subset=["p_nom"])
        .groupby(["country", "carrier"])
        .p_nom.apply(join_exprs)
    )
    minimum = agg_p_nom_minmax["min"].dropna()
    if not minimum.empty:
        minconstraint = define_constraints(
            n, p_nom_per_cc[minimum.index], ">=", minimum, "agg_p_nom", "min"
        )
    maximum = agg_p_nom_minmax["max"].dropna()
    if not maximum.empty:
        maxconstraint = define_constraints(
            n, p_nom_per_cc[maximum.index], "<=", maximum, "agg_p_nom", "max"
        )


def add_EQ_constraints(n, o, scaling=1e-1):
    float_regex = "[0-9]*\.?[0-9]+"
    level = float(re.findall(float_regex, o)[0])
    if o[-1] == "c":
        ggrouper = n.generators.bus.map(n.buses.country)
        lgrouper = n.loads.bus.map(n.buses.country)
        sgrouper = n.storage_units.bus.map(n.buses.country)
    else:
        ggrouper = n.generators.bus
        lgrouper = n.loads.bus
        sgrouper = n.storage_units.bus
    load = (
        n.snapshot_weightings.generators
        @ n.loads_t.p_set.groupby(lgrouper, axis=1).sum()
    )
    inflow = (
        n.snapshot_weightings.stores
        @ n.storage_units_t.inflow.groupby(sgrouper, axis=1).sum()
    )
    inflow = inflow.reindex(load.index).fillna(0.0)
    rhs = scaling * (level * load - inflow)
    lhs_gen = (
        linexpr(
            (n.snapshot_weightings.generators * scaling, get_var(n, "Generator", "p").T)
        )
        .T.groupby(ggrouper, axis=1)
        .apply(join_exprs)
    )
    lhs_spill = (
        linexpr(
            (
                -n.snapshot_weightings.stores * scaling,
                get_var(n, "StorageUnit", "spill").T,
            )
        )
        .T.groupby(sgrouper, axis=1)
        .apply(join_exprs)
    )
    lhs_spill = lhs_spill.reindex(lhs_gen.index).fillna("")
    lhs = lhs_gen + lhs_spill
    define_constraints(n, lhs, ">=", rhs, "equity", "min")


def add_BAU_constraints(n, config):
    ext_c = n.generators.query("p_nom_extendable").carrier.unique()
    mincaps = pd.Series(
        config["electricity"].get("BAU_mincapacities", {key: 0 for key in ext_c})
    )
    lhs = (
        linexpr((1, get_var(n, "Generator", "p_nom")))
        .groupby(n.generators.carrier)
        .apply(join_exprs)
    )
    define_constraints(n, lhs, ">=", mincaps[lhs.index], "Carrier", "bau_mincaps")

    maxcaps = pd.Series(
        config["electricity"].get("BAU_maxcapacities", {key: np.inf for key in ext_c})
    )
    lhs = (
        linexpr((1, get_var(n, "Generator", "p_nom")))
        .groupby(n.generators.carrier)
        .apply(join_exprs)
    )
    define_constraints(n, lhs, "<=", maxcaps[lhs.index], "Carrier", "bau_maxcaps")


def add_SAFE_constraints(n, config):
    peakdemand = (
        1.0 + config["electricity"]["SAFE_reservemargin"]
    ) * n.loads_t.p_set.sum(axis=1).max()
    conv_techs = config["plotting"]["conv_techs"]
    exist_conv_caps = n.generators.query(
        "~p_nom_extendable & carrier in @conv_techs"
    ).p_nom.sum()
    ext_gens_i = n.generators.query("carrier in @conv_techs & p_nom_extendable").index
    lhs = linexpr((1, get_var(n, "Generator", "p_nom")[ext_gens_i])).sum()
    rhs = peakdemand - exist_conv_caps
    define_constraints(n, lhs, ">=", rhs, "Safe", "mintotalcap")


def add_operational_reserve_margin_constraint(n, config):
    reserve_config = config["electricity"]["operational_reserve"]
    EPSILON_LOAD = reserve_config["epsilon_load"]
    EPSILON_VRES = reserve_config["epsilon_vres"]
    CONTINGENCY = reserve_config["contingency"]

    # Reserve Variables
    reserve = get_var(n, "Generator", "r")
    lhs = linexpr((1, reserve)).sum(1)

    # Share of extendable renewable capacities
    ext_i = n.generators.query("p_nom_extendable").index
    vres_i = n.generators_t.p_max_pu.columns
    if not ext_i.empty and not vres_i.empty:
        capacity_factor = n.generators_t.p_max_pu[vres_i.intersection(ext_i)]
        renewable_capacity_variables = get_var(n, "Generator", "p_nom")[
            vres_i.intersection(ext_i)
        ]
        lhs += linexpr(
            (-EPSILON_VRES * capacity_factor, renewable_capacity_variables)
        ).sum(1)

    # Total demand at t
    demand = n.loads_t.p.sum(1)

    # VRES potential of non extendable generators
    capacity_factor = n.generators_t.p_max_pu[vres_i.difference(ext_i)]
    renewable_capacity = n.generators.p_nom[vres_i.difference(ext_i)]
    potential = (capacity_factor * renewable_capacity).sum(1)

    # Right-hand-side
    rhs = EPSILON_LOAD * demand + EPSILON_VRES * potential + CONTINGENCY

    define_constraints(n, lhs, ">=", rhs, "Reserve margin")


def update_capacity_constraint(n):
    gen_i = n.generators.index
    ext_i = n.generators.query("p_nom_extendable").index
    fix_i = n.generators.query("not p_nom_extendable").index

    dispatch = get_var(n, "Generator", "p")
    reserve = get_var(n, "Generator", "r")

    capacity_fixed = n.generators.p_nom[fix_i]

    p_max_pu = get_as_dense(n, "Generator", "p_max_pu")

    lhs = linexpr((1, dispatch), (1, reserve))

    if not ext_i.empty:
        capacity_variable = get_var(n, "Generator", "p_nom")
        lhs += linexpr((-p_max_pu[ext_i], capacity_variable)).reindex(
            columns=gen_i, fill_value=""
        )

    rhs = (p_max_pu[fix_i] * capacity_fixed).reindex(columns=gen_i, fill_value=0)

    define_constraints(n, lhs, "<=", rhs, "Generators", "updated_capacity_constraint")


def add_operational_reserve_margin(n, sns, config):
    """
    Build reserve margin constraints based on the formulation given in
    https://genxproject.github.io/GenX/dev/core/#Reserves.
    """

    define_variables(n, 0, np.inf, "Generator", "r", axes=[sns, n.generators.index])

    add_operational_reserve_margin_constraint(n, config)

    update_capacity_constraint(n)


def add_battery_constraints(n):
    nodes = n.buses.index[n.buses.carrier == "battery"]
    if nodes.empty or ("Link", "p_nom") not in n.variables.index:
        return
    link_p_nom = get_var(n, "Link", "p_nom")
    lhs = linexpr(
        (1, link_p_nom[nodes + " charger"]),
        (
            -n.links.loc[nodes + " discharger", "efficiency"].values,
            link_p_nom[nodes + " discharger"].values,
        ),
    )
    define_constraints(n, lhs, "=", 0, "Link", "charger_ratio")


def add_RES_constraints(n, res_share):
    lgrouper = n.loads.bus.map(n.buses.country)
    ggrouper = n.generators.bus.map(n.buses.country)
    sgrouper = n.storage_units.bus.map(n.buses.country)
    cgrouper = n.links.bus0.map(n.buses.country)

    logger.warning(
        "The add_RES_constraints functionality is still work in progress. "
        "Unexpected results might be incurred, particularly if "
        "temporal clustering is applied or if an unexpected change of technologies "
        "is subject to the obtimisation."
    )

    load = (
        n.snapshot_weightings.generators
        @ n.loads_t.p_set.groupby(lgrouper, axis=1).sum()
    )

    rhs = res_share * load

    res_techs = [
        "solar",
        "onwind",
        "offwind-dc",
        "offwind-ac",
        "battery",
        "hydro",
        "ror",
    ]
    charger = ["H2 electrolysis", "battery charger"]
    discharger = ["H2 fuel cell", "battery discharger"]

    gens_i = n.generators.query("carrier in @res_techs").index
    stores_i = n.storage_units.query("carrier in @res_techs").index
    charger_i = n.links.query("carrier in @charger").index
    discharger_i = n.links.query("carrier in @discharger").index

    # Generators
    lhs_gen = (
        linexpr(
            (n.snapshot_weightings.generators, get_var(n, "Generator", "p")[gens_i].T)
        )
        .T.groupby(ggrouper, axis=1)
        .apply(join_exprs)
    )

    # StorageUnits
    lhs_dispatch = (
        (
            linexpr(
                (
                    n.snapshot_weightings.stores,
                    get_var(n, "StorageUnit", "p_dispatch")[stores_i].T,
                )
            )
            .T.groupby(sgrouper, axis=1)
            .apply(join_exprs)
        )
        .reindex(lhs_gen.index)
        .fillna("")
    )

    lhs_store = (
        (
            linexpr(
                (
                    -n.snapshot_weightings.stores,
                    get_var(n, "StorageUnit", "p_store")[stores_i].T,
                )
            )
            .T.groupby(sgrouper, axis=1)
            .apply(join_exprs)
        )
        .reindex(lhs_gen.index)
        .fillna("")
    )

    # Stores (or their resp. Link components)
    # Note that the variables "p0" and "p1" currently do not exist.
    # Thus, p0 and p1 must be derived from "p" (which exists), taking into account the link efficiency.
    lhs_charge = (
        (
            linexpr(
                (
                    -n.snapshot_weightings.stores,
                    get_var(n, "Link", "p")[charger_i].T,
                )
            )
            .T.groupby(cgrouper, axis=1)
            .apply(join_exprs)
        )
        .reindex(lhs_gen.index)
        .fillna("")
    )

    lhs_discharge = (
        (
            linexpr(
                (
                    n.snapshot_weightings.stores.apply(
                        lambda r: r * n.links.loc[discharger_i].efficiency
                    ),
                    get_var(n, "Link", "p")[discharger_i],
                )
            )
            .groupby(cgrouper, axis=1)
            .apply(join_exprs)
        )
        .reindex(lhs_gen.index)
        .fillna("")
    )

    # signs of resp. terms are coded in the linexpr.
    # todo: for links (lhs_charge and lhs_discharge), account for snapshot weightings
    lhs = lhs_gen + lhs_dispatch + lhs_store + lhs_charge + lhs_discharge

    define_constraints(n, lhs, "=", rhs, "RES share")


def add_land_use_constraint(n):
    if "m" in snakemake.wildcards.clusters:
        _add_land_use_constraint_m(n)
    else:
        _add_land_use_constraint(n)


def _add_land_use_constraint(n):
    # warning: this will miss existing offwind which is not classed AC-DC and has carrier 'offwind'

    for carrier in ["solar", "onwind", "offwind-ac", "offwind-dc"]:
        existing = (
            n.generators.loc[n.generators.carrier == carrier, "p_nom"]
            .groupby(n.generators.bus.map(n.buses.location))
            .sum()
        )
        existing.index += " " + carrier + "-" + snakemake.wildcards.planning_horizons
        n.generators.loc[existing.index, "p_nom_max"] -= existing

    n.generators.p_nom_max.clip(lower=0, inplace=True)


def _add_land_use_constraint_m(n):
    # if generators clustering is lower than network clustering, land_use accounting is at generators clusters

    planning_horizons = snakemake.config["scenario"]["planning_horizons"]
    grouping_years = snakemake.config["existing_capacities"]["grouping_years"]
    current_horizon = snakemake.wildcards.planning_horizons

    for carrier in ["solar", "onwind", "offwind-ac", "offwind-dc"]:
        existing = n.generators.loc[n.generators.carrier == carrier, "p_nom"]
        ind = list(
            set(
                [
                    i.split(sep=" ")[0] + " " + i.split(sep=" ")[1]
                    for i in existing.index
                ]
            )
        )

        previous_years = [
            str(y)
            for y in planning_horizons + grouping_years
            if y < int(snakemake.wildcards.planning_horizons)
        ]

        for p_year in previous_years:
            ind2 = [
                i for i in ind if i + " " + carrier + "-" + p_year in existing.index
            ]
            sel_current = [i + " " + carrier + "-" + current_horizon for i in ind2]
            sel_p_year = [i + " " + carrier + "-" + p_year for i in ind2]
            n.generators.loc[sel_current, "p_nom_max"] -= existing.loc[
                sel_p_year
            ].rename(lambda x: x[:-4] + current_horizon)

    n.generators.p_nom_max.clip(lower=0, inplace=True)


def add_h2_network_cap(n, cap):
    h2_network = n.links.loc[n.links.carrier == "H2 pipeline"]
    if h2_network.index.empty or ("Link", "p_nom") not in n.variables.index:
        return
    h2_network_cap = get_var(n, "Link", "p_nom")
    subset_index = h2_network.index.intersection(h2_network_cap.index)
    lhs = linexpr(
        (h2_network.loc[subset_index, "length"], h2_network_cap[subset_index])
    ).sum()
    # lhs = linexpr((1, h2_network_cap[h2_network.index])).sum()
    rhs = cap * 1000
    define_constraints(n, lhs, "<=", rhs, "h2_network_cap")


def H2_export_yearly_constraint(n):
    res = [
        "csp",
        "rooftop-solar",
        "solar",
        "onwind",
        "onwind2",
        "offwind",
        "offwind2",
        "ror",
    ]
    res_index = n.generators.loc[n.generators.carrier.isin(res)].index

    weightings = pd.DataFrame(
        np.outer(n.snapshot_weightings["generators"], [1.0] * len(res_index)),
        index=n.snapshots,
        columns=res_index,
    )
    res = join_exprs(
        linexpr((weightings, get_var(n, "Generator", "p")[res_index]))
    )  # single line sum

    load_ind = n.loads[n.loads.carrier == "AC"].index.intersection(
        n.loads_t.p_set.columns
    )

    load = (
        n.loads_t.p_set[load_ind].sum(axis=1) * n.snapshot_weightings["generators"]
    ).sum()

    h2_export = n.loads.loc["H2 export load"].p_set * 8760

    lhs = res

    include_country_load = snakemake.config["policy_config"]["yearly"][
        "re_country_load"
    ]

    if include_country_load:
        elec_efficiency = (
            n.links.filter(like="Electrolysis", axis=0).loc[:, "efficiency"].mean()
        )
        rhs = (
            h2_export * (1 / elec_efficiency) + load
        )  # 0.7 is approximation of electrloyzer efficiency # TODO obtain value from network
    else:
        rhs = h2_export * (1 / 0.7)

    con = define_constraints(n, lhs, ">=", rhs, "H2ExportConstraint", "RESproduction")


def monthly_constraints(n, n_ref):
    res_techs = [
        "csp",
        "rooftop-solar",
        "solar",
        "onwind",
        "onwind2",
        "offwind",
        "offwind2",
        "ror",
    ]
    allowed_excess = snakemake.config["policy_config"]["hydrogen"]["allowed_excess"]

    res_index = n.generators.loc[n.generators.carrier.isin(res_techs)].index

    weightings = pd.DataFrame(
        np.outer(n.snapshot_weightings["generators"], [1.0] * len(res_index)),
        index=n.snapshots,
        columns=res_index,
    )

    res = linexpr((weightings, get_var(n, "Generator", "p")[res_index])).sum(
        axis=1
    )  # single line sum
    res = res.groupby(res.index.month).sum()

    electrolysis = get_var(n, "Link", "p")[
        n.links.index[n.links.index.str.contains("H2 Electrolysis")]
    ]
    weightings_electrolysis = pd.DataFrame(
        np.outer(
            n.snapshot_weightings["generators"], [1.0] * len(electrolysis.columns)
        ),
        index=n.snapshots,
        columns=electrolysis.columns,
    )

    elec_input = linexpr((-allowed_excess * weightings_electrolysis, electrolysis)).sum(
        axis=1
    )

    elec_input = elec_input.groupby(elec_input.index.month).sum()

    if snakemake.config["policy_config"]["hydrogen"]["additionality"]:
        res_ref = n_ref.generators_t.p[res_index] * weightings
        res_ref = res_ref.groupby(n_ref.generators_t.p.index.month).sum().sum(axis=1)

        elec_input_ref = (
            n_ref.links_t.p0.loc[
                :, n_ref.links_t.p0.columns.str.contains("H2 Electrolysis")
            ]
            * weightings_electrolysis
        )
        elec_input_ref = (
            -elec_input_ref.groupby(elec_input_ref.index.month).sum().sum(axis=1)
        )

        for i in range(len(res.index)):
            lhs = res.iloc[i] + "\n" + elec_input.iloc[i]
            rhs = res_ref.iloc[i] + elec_input_ref.iloc[i]
            con = define_constraints(
                n, lhs, ">=", rhs, f"RESconstraints_{i}", f"REStarget_{i}"
            )

    else:
        for i in range(len(res.index)):
            lhs = res.iloc[i] + "\n" + elec_input.iloc[i]

            con = define_constraints(
                n, lhs, ">=", 0.0, f"RESconstraints_{i}", f"REStarget_{i}"
            )
    # else:
    #     logger.info("ignoring H2 export constraint as wildcard is set to 0")


def add_chp_constraints(n):
    electric_bool = (
        n.links.index.str.contains("urban central")
        & n.links.index.str.contains("CHP")
        & n.links.index.str.contains("electric")
    )
    heat_bool = (
        n.links.index.str.contains("urban central")
        & n.links.index.str.contains("CHP")
        & n.links.index.str.contains("heat")
    )

    electric = n.links.index[electric_bool]
    heat = n.links.index[heat_bool]

    electric_ext = n.links.index[electric_bool & n.links.p_nom_extendable]
    heat_ext = n.links.index[heat_bool & n.links.p_nom_extendable]

    electric_fix = n.links.index[electric_bool & ~n.links.p_nom_extendable]
    heat_fix = n.links.index[heat_bool & ~n.links.p_nom_extendable]

    link_p = get_var(n, "Link", "p")

    if not electric_ext.empty:
        link_p_nom = get_var(n, "Link", "p_nom")

        # ratio of output heat to electricity set by p_nom_ratio
        lhs = linexpr(
            (
                n.links.loc[electric_ext, "efficiency"]
                * n.links.loc[electric_ext, "p_nom_ratio"],
                link_p_nom[electric_ext],
            ),
            (-n.links.loc[heat_ext, "efficiency"].values, link_p_nom[heat_ext].values),
        )

        define_constraints(n, lhs, "=", 0, "chplink", "fix_p_nom_ratio")

        # top_iso_fuel_line for extendable
        lhs = linexpr(
            (1, link_p[heat_ext]),
            (1, link_p[electric_ext].values),
            (-1, link_p_nom[electric_ext].values),
        )

        define_constraints(n, lhs, "<=", 0, "chplink", "top_iso_fuel_line_ext")

    if not electric_fix.empty:
        # top_iso_fuel_line for fixed
        lhs = linexpr((1, link_p[heat_fix]), (1, link_p[electric_fix].values))

        rhs = n.links.loc[electric_fix, "p_nom"].values

        define_constraints(n, lhs, "<=", rhs, "chplink", "top_iso_fuel_line_fix")

    if not electric.empty:
        # backpressure
        lhs = linexpr(
            (
                n.links.loc[electric, "c_b"].values * n.links.loc[heat, "efficiency"],
                link_p[heat],
            ),
            (-n.links.loc[electric, "efficiency"].values, link_p[electric].values),
        )

        define_constraints(n, lhs, "<=", 0, "chplink", "backpressure")


def add_co2_sequestration_limit(n, sns):
    co2_stores = n.stores.loc[n.stores.carrier == "co2 stored"].index

    if co2_stores.empty or ("Store", "e") not in n.variables.index:
        return

    vars_final_co2_stored = get_var(n, "Store", "e").loc[sns[-1], co2_stores]

    lhs = linexpr((1, vars_final_co2_stored)).sum()
    rhs = (
        n.config["sector"].get("co2_sequestration_potential", 5) * 1e6
    )  # TODO change 200 limit (Europe)

    name = "co2_sequestration_limit"
    define_constraints(
        n, lhs, "<=", rhs, "GlobalConstraint", "mu", axes=pd.Index([name]), spec=name
    )


def set_h2_colors(n):
    blue_h2 = get_var(n, "Link", "p")[
        n.links.index[n.links.index.str.contains("blue H2")]
    ]

    pink_h2 = get_var(n, "Link", "p")[
        n.links.index[n.links.index.str.contains("pink H2")]
    ]

    fuelcell_ind = n.loads[n.loads.carrier == "land transport fuel cell"].index

    other_ind = n.loads[
        (n.loads.carrier == "H2 for industry")
        | (n.loads.carrier == "H2 for shipping")
        | (n.loads.carrier == "H2")
    ].index

    load_fuelcell = (
        n.loads_t.p_set[fuelcell_ind].sum(axis=1) * n.snapshot_weightings["generators"]
    ).sum()

    load_other_h2 = n.loads.loc[other_ind].p_set.sum() * 8760

    load_h2 = load_fuelcell + load_other_h2

    weightings_blue = pd.DataFrame(
        np.outer(n.snapshot_weightings["generators"], [1.0] * len(blue_h2.columns)),
        index=n.snapshots,
        columns=blue_h2.columns,
    )

    weightings_pink = pd.DataFrame(
        np.outer(n.snapshot_weightings["generators"], [1.0] * len(pink_h2.columns)),
        index=n.snapshots,
        columns=pink_h2.columns,
    )

    total_blue = linexpr((weightings_blue, blue_h2)).sum().sum()

    total_pink = linexpr((weightings_pink, pink_h2)).sum().sum()

    rhs_blue = load_h2 * snakemake.config["sector"]["hydrogen"]["blue_share"]
    rhs_pink = load_h2 * snakemake.config["sector"]["hydrogen"]["pink_share"]

    define_constraints(n, total_blue, "=", rhs_blue, "blue_h2_share")

    define_constraints(n, total_pink, "=", rhs_pink, "pink_h2_share")


def add_existing(n):
    if snakemake.wildcards["planning_horizons"] == "2050":
        directory = (
            "results/"
            + "Existing_capacities/"
            + snakemake.config["run"].replace("2050", "2030")
        )
        n_name = (
            snakemake.input.network.split("/")[-1]
            .replace(str(snakemake.config["scenario"]["clusters"][0]), "")
            .replace(str(snakemake.config["costs"]["discountrate"][0]), "")
            .replace("_presec", "")
            .replace(".nc", ".csv")
        )
        df = pd.read_csv(directory + "/electrolyzer_caps_" + n_name, index_col=0)
        existing_electrolyzers = df.p_nom_opt.values

        h2_index = n.links[n.links.carrier == "H2 Electrolysis"].index
        n.links.loc[h2_index, "p_nom_min"] = existing_electrolyzers

        # n_name = snakemake.input.network.split("/")[-1].replace(str(snakemake.config["scenario"]["clusters"][0]), "").\
        #     replace(".nc", ".csv").replace(str(snakemake.config["costs"]["discountrate"][0]), "")
        df = pd.read_csv(directory + "/res_caps_" + n_name, index_col=0)

        for tech in snakemake.config["custom_data"]["renewables"]:
            # df = pd.read_csv(snakemake.config["custom_data"]["existing_renewables"], index_col=0)
            existing_res = df.loc[tech]
            existing_res.index = existing_res.index.str.apply(lambda x: x + tech)
            tech_index = n.generators[n.generators.carrier == tech].index
            n.generators.loc[tech_index, tech] = existing_res


def extra_functionality(n, snapshots):
    """
    Collects supplementary constraints which will be passed to
    ``pypsa.linopf.network_lopf``.

    If you want to enforce additional custom constraints, this is a good location to add them.
    The arguments ``opts`` and ``snakemake.config`` are expected to be attached to the network.
    """
    opts = n.opts
    config = n.config
    if "BAU" in opts and n.generators.p_nom_extendable.any():
        add_BAU_constraints(n, config)
    if "SAFE" in opts and n.generators.p_nom_extendable.any():
        add_SAFE_constraints(n, config)
    if "CCL" in opts and n.generators.p_nom_extendable.any():
        add_CCL_constraints(n, config)
    reserve = config["electricity"].get("operational_reserve", {})
    if reserve.get("activate"):
        add_operational_reserve_margin(n, snapshots, config)
    for o in opts:
        if "RES" in o:
            res_share = float(re.findall("[0-9]*\.?[0-9]+$", o)[0])
            add_RES_constraints(n, res_share)
    for o in opts:
        if "EQ" in o:
            add_EQ_constraints(n, o)
    add_battery_constraints(n)

    if (
        snakemake.config["policy_config"]["hydrogen"]["temporal_matching"]
        == "h2_yearly_matching"
    ):
        if snakemake.config["policy_config"]["hydrogen"]["additionality"] == True:
            logger.info(
                "additionality is currently not supported for yearly constraints, proceeding without additionality"
            )
        logger.info("setting h2 export to yearly greenness constraint")
        H2_export_yearly_constraint(n)

    elif (
        snakemake.config["policy_config"]["hydrogen"]["temporal_matching"]
        == "h2_monthly_matching"
    ):
        if not snakemake.config["policy_config"]["hydrogen"]["is_reference"]:
            logger.info("setting h2 export to monthly greenness constraint")
            monthly_constraints(n, n_ref)
        else:
            logger.info("preparing reference case for additionality constraint")

    elif (
        snakemake.config["policy_config"]["hydrogen"]["temporal_matching"]
        == "no_res_matching"
    ):
        logger.info("no h2 export constraint set")

    else:
        raise ValueError(
            'H2 export constraint is invalid, check config["policy_config"]'
        )

    if snakemake.config["sector"]["hydrogen"]["network"]:
        if snakemake.config["sector"]["hydrogen"]["network_limit"]:
            add_h2_network_cap(
                n, snakemake.config["sector"]["hydrogen"]["network_limit"]
            )

    if snakemake.config["sector"]["hydrogen"]["set_color_shares"]:
        logger.info("setting H2 color mix")
        set_h2_colors(n)

    add_co2_sequestration_limit(n, snapshots)


def solve_network(n, config, solving={}, opts="", **kwargs):
    set_of_options = solving["solver"]["options"]
    cf_solving = solving["options"]

    solver_options = solving["solver_options"][set_of_options] if set_of_options else {}
    solver_name = solving["solver"]["name"]

    track_iterations = cf_solving.get("track_iterations", False)
    min_iterations = cf_solving.get("min_iterations", 4)
    max_iterations = cf_solving.get("max_iterations", 6)

    # add to network for extra_functionality
    n.config = config
    n.opts = opts

    if cf_solving.get("skip_iterations", False):
        network_lopf(
            n,
            solver_name=solver_name,
            solver_options=solver_options,
            extra_functionality=extra_functionality,
            **kwargs,
        )
    else:
        ilopf(
            n,
            solver_name=solver_name,
            solver_options=solver_options,
            track_iterations=track_iterations,
            min_iterations=min_iterations,
            max_iterations=max_iterations,
            extra_functionality=extra_functionality,
            **kwargs,
        )
    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_network",
            simpl="",
            clusters="54",
            ll="copt",
            opts="Co2L-1H",
        )

    configure_logging(snakemake)

    tmpdir = snakemake.params.solving.get("tmpdir")
    if tmpdir is not None:
        Path(tmpdir).mkdir(parents=True, exist_ok=True)
    opts = snakemake.wildcards.opts.split("-")
    solving = snakemake.params.solving

    is_sector_coupled = "sopts" in snakemake.wildcards.keys()

    if is_sector_coupled:
        overrides = override_component_attrs(snakemake.input.overrides)
        n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)
    else:
        n = pypsa.Network(snakemake.input.network)

    if snakemake.params.augmented_line_connection.get("add_to_snakefile"):
        n.lines.loc[n.lines.index.str.contains("new"), "s_nom_min"] = (
            snakemake.params.augmented_line_connection.get("min_expansion")
        )

    if (
        snakemake.config["custom_data"]["add_existing"]
        and snakemake.wildcards.planning_horizons == "2050"
        and is_sector_coupled
    ):
        add_existing(n)

    if (
        snakemake.config["policy_config"]["hydrogen"]["additionality"]
        and not snakemake.config["policy_config"]["hydrogen"]["is_reference"]
        and snakemake.config["policy_config"]["hydrogen"]["temporal_matching"]
        != "no_res_matching"
        and is_sector_coupled
    ):
        n_ref_path = snakemake.config["policy_config"]["hydrogen"]["path_to_ref"]
        n_ref = pypsa.Network(n_ref_path)
    else:
        n_ref = None

    n = prepare_network(n, solving["options"])

    n = solve_network(
        n,
        config=snakemake.config,
        solving=solving,
        opts=opts,
        solver_dir=tmpdir,
        solver_logfile=snakemake.log.solver,
    )
    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output[0])
    logger.info(f"Objective function: {n.objective}")
    logger.info(f"Objective constant: {n.objective_constant}")
