# -*- coding: utf-8 -*-

import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pypsa
import pytz
import xarray as xr
from helpers import (
    create_dummy_data,
    create_network_topology,
    cycling_shift,
    locate_bus,
    mock_snakemake,
    override_component_attrs,
    prepare_costs,
    three_2_two_digits_country,
    two_2_three_digits_country,
)

if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        # from helper import mock_snakemake #TODO remove func from here to helper script
        snakemake = mock_snakemake(
            "prepare_sector_network",
            simpl="",
            clusters="53",
            ll="c1",
            opts="",
            planning_horizons="2030",
            sopts="Co2L-720H",
        )

    options = ["solar", "offwind", "onwind", "csp"]
    options = ["solar"]  # , 'offwind', 'onwind', 'csp']

    n = pypsa.Network(
        "../results/v.0.0.1_find_bug6/prenetworks/elec_s_285_ec_lc1.0_Co2L_730H_2030.nc"
    )

    # for option in options:

    #     if option in n.generators.carrier:
    #         gen_ind = n.generators[n.generators.carrier=='solar'].index
    #         n.generators[gen_ind]= custom_res['']
    #     else
