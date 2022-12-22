# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2022 PyPSA-Earth Authors
#
# coding: utf-8
"""
Execute a scenario optimization

Run iteratively the workflow under different conditions
and store the results in specific folders
"""
import collections.abc
import copy
import os
import shutil

from _helpers import sets_path_to_root
from build_test_configs import create_test_config

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("run_scenario", scenario="NG")

    sets_path_to_root("pypsa-earth")

    scenario = snakemake.wildcards["scenario"]
    base_config = snakemake.config.get("base_config", "./config.default.yaml")

    # create scenario config
    create_test_config(base_config, f"configs/scenarios/{scenario}.yaml", "config.yaml")

    val = os.system("snakemake -j all solve_all_networks --forceall")

    for f in ["resources", "networks", "results"]:
        shutil.copytree(f, f"scenarios/{scenario}/{f}")

    shutil.copy("config.yaml", f"scenarios/{scenario}/config.yaml")