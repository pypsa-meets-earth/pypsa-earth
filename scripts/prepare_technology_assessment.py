# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2022 The PyPSA-Earth Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later
# coding: utf-8
"""

"""
import logging
import os
import numpy as np
import pandas as pd
import pypsa
from _helpers import configure_logging

logger = logging.getLogger(__name__)


def single_best_in_worst_list(worst_list, best_list):
    """
    Return list with single best value per list of worst values
    
    Input:
    ------
    worst_list: 1D array
    best_list: 1D array
    
    Output:
    -------
    new_list: 2D array
    
    Example:
    --------
    >>> single_best_in_worst_list([1, 1, 1], [4, 4, 4])
    [[4, 1, 1], [1, 4, 1], [1, 1, 4]]
    """
    new_list = []
    for i in range(len(worst_list)):
        l = worst_list.copy()
        l[i] = best_list[i]
        new_list.append(l)
    
    return new_list


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "prepare_technology_assessment",
            simpl="",
            clusters="6",
            ll="copt",
            opts="Co2L-4H",
            unc="f0",
        )
    configure_logging(snakemake)
    CONFIG = snakemake.config
    n = pypsa.Network(snakemake.input[0])

    ### SCENARIO INPUTS
    ###
    PYPSA_FEATURES = {
        k: v for k, v in CONFIG["technology_assessment"]["pypsa_standard"].items() if v
    }  # removes key value pairs with empty value e.g. []
    L_BOUNDS = [item[0] for item in PYPSA_FEATURES.values()]
    U_BOUNDS = [item[1] for item in PYPSA_FEATURES.values()]
    STRATEGY = CONFIG["technology_assessment"]["options"].get("assessment_strategy")
    

    if STRATEGY=="any-chance":
        carrier_no = len(n.stores.carrier.unique())
        worst_list = L_BOUNDS * carrier_no
        best_list = U_BOUNDS * carrier_no
        scenarios = single_best_in_worst_list(worst_list, best_list)

    ### NETWORK CHANGES
    ###
    unc_wildcards = snakemake.wildcards[-1]
    i = int(unc_wildcards[1:])
    for k, _ in PYPSA_FEATURES.items():
        type = k.split(".")[0] # "stores", "generators", ...
        feature = k.split(".")[1] # "capital_cost", "efficiency", ...

    if type == "stores":
        carrier_list = n.stores.carrier.unique()
        j = 0
        for c in carrier_list:
            n.stores.loc[n.stores.carrier==c, feature] *= scenarios[i][j]
            n.links.loc[n.links.carrier.str.contains("H2"), feature] *= scenarios[i][j] 
            logger.info(f"Scaled {feature} for carrier={c} of store and links by factor {scenarios[i][j]} in the {i} scenario")
            j += 1

    ### EXPORT AND METADATA
    #
    scenario_dict = (
        pd.DataFrame(scenarios).rename_axis("iteration").add_suffix("_carrier")
    ).to_dict()
    n.meta.update(scenario_dict)
    n.export_to_netcdf(snakemake.output[0])
