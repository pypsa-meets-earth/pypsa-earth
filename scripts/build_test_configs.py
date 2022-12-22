# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2022 PyPSA-Earth Authors
#
# coding: utf-8
"""
Write option files (configs) for the Continous Integration tests

The config.tutorial.yaml has all options.
The test/* config files have only key/value strings that are different from the tutorial config.
The below scripts 'updates' the test configs and adds all options of the tutorial config.
"""
import collections.abc
import copy
import os
from pathlib import Path

from ruamel.yaml import YAML


def update(d, u):
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def create_test_config(base_path, changes_path, output_path):
    # Input paths
    fp_baseconfig = Path(Path.cwd(), base_path)
    fp_update_file = Path(Path.cwd(), changes_path)

    # Load yaml files
    yaml = YAML()
    with open(fp_baseconfig) as fp:
        base_config = yaml.load(fp)

    # Load update yaml
    with open(fp_update_file) as fp:
        update_config = yaml.load(fp)
    # create updated yaml
    test_config = update(copy.deepcopy(base_config), update_config)
    # Output path
    fp = Path(Path.cwd(), output_path)
    # Save file
    yaml.dump(test_config, fp)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_test_configs")

    # Input paths
    fp_baseconfig = snakemake.input.base_config
    fp_update_file_list = snakemake.input.update_file_list
    fp_output_file_list = snakemake.output.tmp_test_configs

    for (finput, foutput) in zip(fp_update_file_list, fp_output_file_list):
        create_test_config(fp_baseconfig, finput, foutput)

    # Manual output in terminal
    # import sys
    # yaml.dump(data, sys.stdout)
