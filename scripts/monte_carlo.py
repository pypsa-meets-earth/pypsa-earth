# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2022 The PyPSA-Africa Authors
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

import chaospy
import numpy as np
import pandas as pd
import pypsa
from _helpers import configure_logging
from pyDOE2 import lhs
from scipy.stats import qmc
from solve_network import *

logger = logging.getLogger(__name__)


def monte_carlo_sampling_pydoe2(
    N_FEATURES,
    SAMPLES,
    criterion=None,
    iteration=None,
    random_state=42,
    correlation_matrix=None,
):
    """
    Creates Latin Hypercube Sample (LHS) implementation from PyDOE2 with various options. Additionally all "corners" are simulated.

    Adapted from Disspaset: https://github.com/energy-modelling-toolkit/Dispa-SET/blob/master/scripts/build_and_run_hypercube.py
    Documentation on PyDOE2: https://github.com/clicumu/pyDOE2 (fixes latin_cube errors)
    """

    # Generate a Nfeatures-dimensional latin hypercube varying between 0 and 1:
    lh = lhs(
        N_FEATURES,
        samples=SAMPLES,
        criterion=criterion,
        iterations=iteration,
        random_state=random_state,
        correlation_matrix=correlation_matrix,
    )
    discrepancy = qmc.discrepancy(lh)
    logger.info("Hypercube discrepancy is:", discrepancy)

    return lh


def monte_carlo_sampling_chaospy(N_FEATURES, SAMPLES, rule="latin_hypercube", seed=42):
    """
    Creates Latin Hypercube Sample (LHS) implementation from chaospy.

    Documentation on Chaospy: https://github.com/clicumu/pyDOE2 (fixes latin_cube errors)
    Documentation on Chaospy latin-hyper cube (quasi-Monte Carlo method): https://chaospy.readthedocs.io/en/master/user_guide/fundamentals/quasi_random_samples.html#Quasi-random-samples

    """
    # Generate a Nfeatures-dimensional latin hypercube varying between 0 and 1:
    N_FEATURES = "chaospy.Uniform(0, 1), " * N_FEATURES
    uniform_cube = eval(
        f"chaospy.J({N_FEATURES})"
    )  # writes Nfeatures times the chaospy.uniform... command)
    lh = uniform_cube.sample(SAMPLES, rule=rule, seed=seed).T
    discrepancy = qmc.discrepancy(lh)
    logger.info("Hypercube discrepancy is:", discrepancy)

    return lh


def monte_carlo_sampling_scipy(
    N_FEATURES, SAMPLES, centered=False, strength=2, optimization=None, seed=42
):
    """
    Creates Latin Hypercube Sample (LHS) implementation from SciPy with various options:
    - Center the point within the multi-dimensional grid, centered=True
    - optimization scheme, optimization="random-cd"
    - strength=1, classical LHS
    - strength=2, performant orthogonal LHS, requires the sample to be a prime or square of a prime e.g. sqr(121)=11

    Options could be combined to produce an optimized centered orthogonal array
    based LHS. After optimization, the result would not be guaranteed to be of strength 2.

    Documentation for Quasi-Monte Carlo approaches: https://docs.scipy.org/doc/scipy/reference/stats.qmc.html
    Documentation for Latin Hypercube: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.qmc.LatinHypercube.html#scipy.stats.qmc.LatinHypercube
    Orthogonal LHS is better than basic LHS: https://github.com/scipy/scipy/pull/14546/files, https://en.wikipedia.org/wiki/Latin_hypercube_sampling
    """
    sampler = qmc.LatinHypercube(
        d=N_FEATURES,
        centered=centered,
        strength=strength,
        optimization=optimization,
        seed=seed,
    )
    lh = sampler.random(n=SAMPLES)
    discrepancy = qmc.discrepancy(lh)
    logger.info("Hypercube discrepancy is:", discrepancy)

    return lh


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "monte_carlo",
            network="elec",
            simpl="",
            clusters="20",
            ll="v0.3",
            opts="Co2L-24H",
        )
    configure_logging(snakemake)
    config = snakemake.config

    ### SCENARIO INPUTS
    ###
    MONTE_CARLO_PYPSA_FEATURES = {
        k: v for k, v in config["monte_carlo"]["pypsa_standard"].items() if v
    }  # removes key value pairs with empty value e.g. []
    MONTE_CARLO_OPTIONS = config["monte_carlo"]["options"]
    L_BOUNDS = [item[0] for item in MONTE_CARLO_PYPSA_FEATURES.values()]
    U_BOUNDS = [item[1] for item in MONTE_CARLO_PYPSA_FEATURES.values()]
    N_FEATURES = len(
        MONTE_CARLO_PYPSA_FEATURES
    )  # only counts features when specified in config
    SAMPLES = MONTE_CARLO_OPTIONS.get(
        "samples"
    )  # What is the optimal sampling? Probably depend on amount of features
    SAMPLING_STRATEGY = MONTE_CARLO_OPTIONS.get("sampling_strategy")

    ### SCENARIO CREATION / SAMPLING STRATEGY
    ###
    if SAMPLING_STRATEGY == "pydoe2":
        lh = monte_carlo_sampling_pydoe2(
            N_FEATURES,
            SAMPLES,
            criterion=None,
            iteration=None,
            random_state=42,
            correlation_matrix=None,
        )
    if SAMPLING_STRATEGY == "scipy":
        lh = monte_carlo_sampling_scipy(
            N_FEATURES, SAMPLES, centered=False, strength=2, optimization=None, seed=42
        )
    if SAMPLING_STRATEGY == "chaospy":
        lh = monte_carlo_sampling_chaospy(
            N_FEATURES, SAMPLES, rule="latin_hypercube", seed=42
        )
    lh_scaled = qmc.scale(lh, L_BOUNDS, U_BOUNDS)

    ### SCENARIO ITERATION
    ###
    ### solve_network.py preparation
    tmpdir = config["solving"].get("tmpdir")
    if tmpdir is not None:
        Path(tmpdir).mkdir(parents=True, exist_ok=True)
    opts = snakemake.wildcards.opts.split("-")
    solve_opts = config["solving"]["options"]

    fn = getattr(snakemake.log, "memory", None)
    with memory_logger(filename=fn, interval=30.0) as mem:
        n = pypsa.Network(snakemake.input[0])
        n = prepare_network(n, solve_opts)

    ### monte-carlo part
    Nruns = lh.shape[0]
    for i in range(Nruns):
        nm = n.copy()
        j = 0
        for k, v in MONTE_CARLO_PYPSA_FEATURES.items():
            # this loop sets in one scenario each "i" feature assumption
            # i is the config input key "loads_t.p_set"
            # v is the lower and upper bound [0.8,1.3]
            # j is the sample interation number
            # Example: n.loads_t.p_set = network.loads_t.p_set = .loads_t.p_set * lh[0,0] * (1.3-0.8)
            exec(f"n.{k} = network.{k} * {lh_scaled[i,j]}")
            logger.info(f"Scaled n.{k} by factor {lh_scaled[i,j]} in the {i} scenario")
            j = j + 1

        # run optimization
        nm = solve_network(
            n,
            config=config,
            opts=opts,
            solver_dir=tmpdir,
            solver_logfile=snakemake.log.solver,
        )

        # create new object 'monte_carlo' to save settings for each run
        nm.monte_carlo = (
            pd.DataFrame(lh_scaled).rename_axis("Nruns").add_suffix("_feature")
        )
        nm.export_to_netcdf(
            snakemake.output[0]
        )  #### TODO... Create outputs & parallize

        logger.info("Maximum memory usage: {}".format(mem.mem_usage))
