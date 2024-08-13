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
The optimization is based on the :func:`network.optimize` function.
Additionally, some extra constraints specified in :mod:`prepare_network` and :mod:`solve_network` are added.

Solving the network in multiple iterations is motivated through the dependence of transmission line capacities and impedances on values of corresponding flows.
As lines are expanded their electrical parameters change, which renders the optimisation bilinear even if the power flow
equations are linearized.
To retain the computational advantage of continuous linear programming, a sequential linear programming technique
is used, where in between iterations the line impedances are updated.
Details (and errors introduced through this heuristic) are discussed in the paper

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
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
from _helpers import configure_logging, create_logger
from linopy import merge
from pypsa.descriptors import get_switchable_as_dense as get_as_dense
from pypsa.optimization.abstract import optimize_transmission_expansion_iteratively
from pypsa.optimization.optimize import optimize
from vresutils.benchmark import memory_logger

logger = create_logger(__name__)


def prepare_network(n, solve_opts, config):
    if "clip_p_max_pu" in solve_opts:
        for df in (n.generators_t.p_max_pu, n.storage_units_t.inflow):
            df.where(df > solve_opts["clip_p_max_pu"], other=0.0, inplace=True)

    load_shedding = solve_opts.get("load_shedding")
    if load_shedding:
        n.add("Carrier", "Load")
        buses_i = n.buses.query("carrier == 'AC'").index
        if not np.isscalar(load_shedding):
            load_shedding = 8e3  # Eur/kWh
        # intersect between macroeconomic and surveybased
        # willingness to pay
        # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full)
        # 1e2 is practical relevant, 8e3 good for debugging
        n.madd(
            "Generator",
            buses_i,
            " load",
            bus=buses_i,
            carrier="load",
            sign=1e-3,  # Adjust sign to measure p and p_nom in kW instead of MW
            marginal_cost=load_shedding,
            p_nom=1e9,  # kW
        )

    if solve_opts.get("noisy_costs"):
        for t in n.iterate_components(n.one_port_components):
            # TODO: uncomment out to and test noisy_cost (makes solution unique)
            # if 'capital_cost' in t.df:
            #    t.df['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.df)) - 0.5)
            if "marginal_cost" in t.df:
                t.df["marginal_cost"] += 1e-2 + 2e-3 * (
                    np.random.random(len(t.df)) - 0.5
                )

        for t in n.iterate_components(["Line", "Link"]):
            t.df["capital_cost"] += (
                1e-1 + 2e-2 * (np.random.random(len(t.df)) - 0.5)
            ) * t.df["length"]

    if solve_opts.get("nhours"):
        nhours = solve_opts["nhours"]
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760.0 / nhours

    return n


def add_CCL_constraints(n, config):
    """
    Add CCL (country & carrier limit) constraint to the network.

    Add minimum and maximum levels of generator nominal capacity per carrier
    for individual countries. Opts and path for agg_p_nom_minmax.csv must be defined
    in config.yaml. Default file is available at data/agg_p_nom_minmax.csv.

    Parameters
    ----------
    n : pypsa.Network
    config : dict

    Example
    -------
    scenario:
        opts: [Co2L-CCL-24H]
    electricity:
        agg_p_nom_limits: data/agg_p_nom_minmax.csv
    """
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
    capacity_variable = n.model["Generator-p_nom"]

    lhs = []
    ext_carriers = n.generators.query("p_nom_extendable").carrier.unique()
    for c in ext_carriers:
        ext_carrier = n.generators.query("p_nom_extendable and carrier == @c")
        country_grouper = (
            ext_carrier.bus.map(n.buses.country)
            .rename_axis("Generator-ext")
            .rename("country")
        )
        ext_carrier_per_country = capacity_variable.loc[
            country_grouper.index
        ].groupby_sum(country_grouper)
        lhs.append(ext_carrier_per_country)
    lhs = merge(lhs, dim=pd.Index(ext_carriers, name="carrier"))

    min_matrix = agg_p_nom_minmax["min"].to_xarray().unstack().reindex_like(lhs)
    max_matrix = agg_p_nom_minmax["max"].to_xarray().unstack().reindex_like(lhs)

    n.model.add_constraints(
        lhs >= min_matrix, name="agg_p_nom_min", mask=min_matrix.notnull()
    )
    n.model.add_constraints(
        lhs <= max_matrix, name="agg_p_nom_max", mask=max_matrix.notnull()
    )


def add_EQ_constraints(n, o, scaling=1e-1):
    """
    Add equity constraints to the network.

    Currently this is only implemented for the electricity sector only.

    Opts must be specified in the config.yaml.

    Parameters
    ----------
    n : pypsa.Network
    o : str

    Example
    -------
    scenario:
        opts: [Co2L-EQ0.7-24h]

    Require each country or node to on average produce a minimal share
    of its total electricity consumption itself. Example: EQ0.7c demands each country
    to produce on average at least 70% of its consumption; EQ0.7 demands
    each node to produce on average at least 70% of its consumption.
    """
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
    dispatch_variable = n.model["Generator-p"]
    lhs_gen = (
        (dispatch_variable * (n.snapshot_weightings.generators * scaling))
        .groupby(ggrouper.to_xarray())
        .sum()
        .sum("snapshot")
    )
    # TODO: double check that this is really needed, why do have to subtract the spillage
    if not n.storage_units_t.inflow.empty:
        spillage_variable = n.model["StorageUnit-spill"]
        lhs_spill = (
            (spillage_variable * (-n.snapshot_weightings.stores * scaling))
            .groupby_sum(sgrouper)
            .groupby(sgrouper.to_xarray())
            .sum()
            .sum("snapshot")
        )
        lhs = lhs_gen + lhs_spill
    else:
        lhs = lhs_gen
    n.model.add_constraints(lhs >= rhs, name="equity_min")


def add_BAU_constraints(n, config):
    """
    Add a per-carrier minimal overall capacity.

    BAU_mincapacities and opts must be adjusted in the config.yaml.

    Parameters
    ----------
    n : pypsa.Network
    config : dict

    Example
    -------
    scenario:
        opts: [Co2L-BAU-24h]
    electricity:
        BAU_mincapacities:
            solar: 0
            onwind: 0
            OCGT: 100000
            offwind-ac: 0
            offwind-dc: 0
    Which sets minimum expansion across all nodes e.g. in Europe to 100GW.
    OCGT bus 1 + OCGT bus 2 + ... > 100000
    """
    mincaps = pd.Series(config["electricity"]["BAU_mincapacities"])
    p_nom = n.model["Generator-p_nom"]
    ext_i = n.generators.query("p_nom_extendable")
    ext_carrier_i = xr.DataArray(ext_i.carrier.rename_axis("Generator-ext"))
    lhs = p_nom.groupby(ext_carrier_i).sum()
    rhs = mincaps[lhs.indexes["carrier"]].rename_axis("carrier")
    n.model.add_constraints(lhs >= rhs, name="bau_mincaps")


def add_SAFE_constraints(n, config):
    """
    Add a capacity reserve margin of a certain fraction above the peak demand.
    Renewable generators and storage do not contribute. Ignores network.

    Parameters
    ----------
        n : pypsa.Network
        config : dict

    Example
    -------
    config.yaml requires to specify opts:

    scenario:
        opts: [Co2L-SAFE-24h]
    electricity:
        SAFE_reservemargin: 0.1
    Which sets a reserve margin of 10% above the peak demand.
    """
    peakdemand = n.loads_t.p_set.sum(axis=1).max()
    margin = 1.0 + config["electricity"]["SAFE_reservemargin"]
    reserve_margin = peakdemand * margin
    conventional_carriers = config["electricity"]["conventional_carriers"]
    ext_gens_i = n.generators.query(
        "carrier in @conventional_carriers & p_nom_extendable"
    ).index
    capacity_variable = n.model["Generator-p_nom"]
    p_nom = n.model["Generator-p_nom"].loc[ext_gens_i]
    lhs = p_nom.sum()
    exist_conv_caps = n.generators.query(
        "~p_nom_extendable & carrier in @conventional_carriers"
    ).p_nom.sum()
    rhs = reserve_margin - exist_conv_caps
    n.model.add_constraints(lhs >= rhs, name="safe_mintotalcap")


def add_operational_reserve_margin_constraint(n, sns, config):
    """
    Build reserve margin constraints based on the formulation
    as suggested in GenX
    https://energy.mit.edu/wp-content/uploads/2017/10/Enhanced-Decision-Support-for-a-Changing-Electricity-Landscape.pdf
    It implies that the reserve margin also accounts for optimal
    dispatch of distributed energy resources (DERs) and demand response
    which is a novel feature of GenX.
    """
    reserve_config = config["electricity"]["operational_reserve"]
    EPSILON_LOAD = reserve_config["epsilon_load"]
    EPSILON_VRES = reserve_config["epsilon_vres"]
    CONTINGENCY = reserve_config["contingency"]

    # Reserve Variables
    n.model.add_variables(
        0, np.inf, coords=[sns, n.generators.index], name="Generator-r"
    )
    reserve = n.model["Generator-r"]
    lhs = reserve.sum("Generator")

    # Share of extendable renewable capacities
    ext_i = n.generators.query("p_nom_extendable").index
    vres_i = n.generators_t.p_max_pu.columns
    if not ext_i.empty and not vres_i.empty:
        capacity_factor = n.generators_t.p_max_pu[vres_i.intersection(ext_i)]
        renewable_capacity_variables = (
            n.model["Generator-p_nom"]
            .loc[vres_i.intersection(ext_i)]
            .rename({"Generator-ext": "Generator"})
        )
        lhs = merge(
            lhs,
            (renewable_capacity_variables * (-EPSILON_VRES * capacity_factor)).sum(
                ["Generator"]
            ),
        )

    # Total demand per t
    demand = get_as_dense(n, "Load", "p_set").sum(axis=1)

    # VRES potential of non extendable generators
    capacity_factor = n.generators_t.p_max_pu[vres_i.difference(ext_i)]
    renewable_capacity = n.generators.p_nom[vres_i.difference(ext_i)]
    potential = (capacity_factor * renewable_capacity).sum(axis=1)

    # Right-hand-side
    rhs = EPSILON_LOAD * demand + EPSILON_VRES * potential + CONTINGENCY

    n.model.add_constraints(lhs >= rhs, name="reserve_margin")


def update_capacity_constraint(n):
    gen_i = n.generators.index
    ext_i = n.generators.query("p_nom_extendable").index
    fix_i = n.generators.query("not p_nom_extendable").index

    dispatch = n.model["Generator-p"]
    reserve = n.model["Generator-r"]
    capacity_fixed = n.generators.p_nom[fix_i]

    p_max_pu = get_as_dense(n, "Generator", "p_max_pu")

    lhs = merge(
        dispatch * 1,
        reserve * 1,
    )

    if not ext_i.empty:
        capacity_variable = n.model["Generator-p_nom"]
        lhs = dispatch + reserve - capacity_variable * xr.DataArray(p_max_pu[ext_i])

    rhs = (p_max_pu[fix_i] * capacity_fixed).reindex(columns=gen_i, fill_value=0)
    n.model.add_constraints(
        lhs <= rhs, name="gen_updated_capacity_constraint", mask=rhs.notnull()
    )


def add_operational_reserve_margin(n, sns, config):
    """
    Build reserve margin constraints based on the formulation given in
    https://genxproject.github.io/GenX/dev/core/#Reserves.

    Parameters
    ----------
        n : pypsa.Network
        sns: pd.DatetimeIndex
        config : dict

    Example:
    --------
    config.yaml requires to specify operational_reserve:
    operational_reserve: # like https://genxproject.github.io/GenX/dev/core/#Reserves
        activate: true
        epsilon_load: 0.02 # percentage of load at each snapshot
        epsilon_vres: 0.02 # percentage of VRES at each snapshot
        contingency: 400000 # MW
    """
    add_operational_reserve_margin_constraint(n, sns, config)
    update_capacity_constraint(n)


def add_battery_constraints(n):
    """
    Add constraint ensuring that charger = discharger, i.e.
    1 * charger_size - efficiency * discharger_size = 0
    """
    nodes = n.buses.index[n.buses.carrier == "battery"]
    # if nodes.empty or ("Link", "p_nom") not in n.variables.index:
    if nodes.empty:
        return
    vars_link = n.model["Link-p_nom"]
    eff = n.links.loc[nodes + " discharger", "efficiency"]
    lhs = merge(
        vars_link.sel({"Link-ext": nodes + " charger"}) * 1,
        # for some reasons, eff is one element longer as compared with vars_link
        vars_link.sel({"Link-ext": nodes + " discharger"}) * -eff[0],
    )
    n.model.add_constraints(lhs == 0, name="link_charger_ratio")


def add_RES_constraints(n, res_share, config):
    lgrouper = n.loads.bus.map(n.buses.country)
    # TODO drop load
    ggrouper = n.generators.bus.map(n.buses.country)
    sgrouper = n.storage_units.bus.map(n.buses.country)
    cgrouper = n.links.bus0.map(n.buses.country)

    logger.warning(
        "The add_RES_constraints() is still work in progress. "
        "Unexpected results might be incurred, particularly if "
        "temporal clustering is applied or if an unexpected change of technologies "
        "is subject to the optimisation."
    )

    load = (
        n.snapshot_weightings.generators
        @ n.loads_t.p_set.groupby(lgrouper, axis=1).sum()
    )

    rhs = res_share * load

    renew_techs = config["electricity"]["renewable_carriers"]
    charger = ["H2 electrolysis", "battery charger"]
    discharger = ["H2 fuel cell", "battery discharger"]

    gens_i = n.generators.query("carrier in @renew_techs").index

    stores_i = n.storage_units.query("carrier in @renew_techs").index
    charger_i = n.links.query("carrier in @charger").index
    discharger_i = n.links.query("carrier in @discharger").index

    stores_t_weights = n.snapshot_weightings.stores

    # Generators
    # TODO restore grouping by countries un-commenting calls of groupby()
    lhs_gen = (
        (n.model["Generator-p"].loc[:, gens_i] * n.snapshot_weightings.generators)
        # .groupby(ggrouper.to_xarray())
        .sum()
    )

    store_disp_expr = (
        n.model["StorageUnit-p_dispatch"].loc[:, stores_i] * stores_t_weights
    )
    store_expr = n.model["StorageUnit-p_store"].loc[:, stores_i] * stores_t_weights
    charge_expr = n.model["Link-p"].loc[:, charger_i] * stores_t_weights.apply(
        lambda r: r * n.links.loc[charger_i].efficiency
    )
    discharge_expr = n.model["Link-p"].loc[:, discharger_i] * stores_t_weights.apply(
        lambda r: r * n.links.loc[discharger_i].efficiency
    )

    lhs_dispatch = (
        store_disp_expr
        # .groupby(sgrouper)
        .sum()
    )
    lhs_store = (
        store_expr
        # .groupby(sgrouper)
        .sum()
    )
    lhs_charge = (
        charge_expr
        # .groupby(cgrouper)
        .sum()
    )
    lhs_discharge = (
        discharge_expr
        # .groupby(cgrouper)
        .sum()
    )

    lhs = lhs_gen + lhs_dispatch - lhs_store - lhs_charge + lhs_discharge

    n.model.add_constraints(lhs == rhs, name="res_share")


def extra_functionality(n, snapshots):
    """
    Collects supplementary constraints which will be passed to
    ``pypsa.optimization.optimize``.

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
            add_RES_constraints(n, res_share, config)
    for o in opts:
        if "EQ" in o:
            add_EQ_constraints(n, o)
    add_battery_constraints(n)


def solve_network(n, config, solving, **kwargs):
    set_of_options = solving["solver"]["options"]
    cf_solving = solving["options"]

    kwargs["solver_options"] = (
        solving["solver_options"][set_of_options] if set_of_options else {}
    )
    kwargs["solver_name"] = solving["solver"]["name"]

    skip_iterations = cf_solving.get("skip_iterations", False)
    if not n.lines.s_nom_extendable.any():
        skip_iterations = True
        logger.info("No expandable lines found. Skipping iterative solving.")

    # add to network for extra_functionality
    n.config = config

    if skip_iterations:
        status, condition = n.optimize(**kwargs)
    else:
        kwargs["track_iterations"] = (cf_solving.get("track_iterations", False),)
        kwargs["min_iterations"] = (cf_solving.get("min_iterations", 4),)
        kwargs["max_iterations"] = (cf_solving.get("max_iterations", 6),)
        status, condition = n.optimize.optimize_transmission_expansion_iteratively(
            **kwargs
        )

    if status != "ok":  # and not rolling_horizon:
        logger.warning(
            f"Solving status '{status}' with termination condition '{condition}'"
        )
    if "infeasible" in condition:
        labels = n.model.compute_infeasibilities()
        logger.info(f"Labels:\n{labels}")
        n.model.print_infeasibilities()
        raise RuntimeError("Solving status 'infeasible'")

    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "solve_network",
            simpl="",
            clusters="54",
            ll="copt",
            opts="Co2L-1H",
        )
    configure_logging(snakemake)

    opts = snakemake.wildcards.opts.split("-")
    solve_opts = snakemake.params.solving["options"]

    n = pypsa.Network(snakemake.input[0])
    # needed to get `n.model` property
    n.optimize.create_model()
    # TODO Double-check handling the augmented case
    # if snakemake.params.augmented_line_connection.get("add_to_snakefile"):
    #    n.lines.loc[
    #        n.lines.index.str.contains("new"), "s_nom_min"
    #    ] = snakemake.params.augmented_line_connection.get("min_expansion")
    n = prepare_network(n, solve_opts, config=solve_opts)
    n = solve_network(
        n,
        config=snakemake.config,
        solving=snakemake.params.solving,
        log_fn=snakemake.log.solver,
    )
    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output[0])
    logger.info(f"Objective function: {n.objective}")
    logger.info(f"Objective constant: {n.objective_constant}")
