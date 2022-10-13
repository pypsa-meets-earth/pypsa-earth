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


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_test_configs")

    # Input paths
    fp_baseconfig = Path(Path.cwd(), snakemake.input.base_config)
    fp_update_file_list = [Path(Path.cwd(), i) for i in snakemake.input.update_file_list]

    # Load yaml files
    yaml = YAML()
    with open(fp_baseconfig) as fp:
        baseconfig = yaml.load(fp)
    
    with open(fp_update_file_list[0]) as fp:
        base_update = yaml.load(fp_update_file_list[0])
        base_test_config = update(copy.deepcopy(baseconfig), base_update)
        fp = Path(Path.cwd(), snakemake.output.tmp_test_configs[0])
        yaml.dump(base_test_config, fp)

    for c in fp_update_file_list[1:]:
        # Load update yaml
        with open(c) as fp:
            update_config = yaml.load(fp)
        # create updated yaml
        test_config = update(copy.deepcopy(base_test_config), update_config)
        # Output path
        list_no = fp_update_file_list.index(c)
        fp = Path(Path.cwd(), snakemake.output[list_no])
        # Save file
        yaml.dump(test_config, fp)


    # Manual output in terminal
    # import sys
    # yaml.dump(data, sys.stdout)
