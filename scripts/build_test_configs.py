# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Write option files (configs) for the Continuous Integration tests.

The config.tutorial.yaml has all options. The test/* config files have
only key/value strings that are different from the tutorial config. The
below scripts 'updates' the test configs and adds all options of the
tutorial config.
"""
import collections.abc
import copy
import os
from pathlib import Path

from ruamel.yaml import YAML


def update(base_dictionary, diff_dictionary):
    """
    A function to recursivelly update items in a dictionary which
    intended usage is unpdating the base config according to the diff
    config
    """
    for key_to_apply, value_to_apply in diff_dictionary.items():
        if isinstance(value_to_apply, collections.abc.Mapping) and isinstance(
            base_dictionary.get(key_to_apply), collections.abc.Mapping
        ):
            base_keys = set(base_dictionary[key_to_apply])
            diff_keys = set(value_to_apply)

            if base_keys & diff_keys:
                base_dictionary[key_to_apply] = update(
                    base_dictionary[key_to_apply], value_to_apply
                )
            else:
                base_dictionary[key_to_apply] = copy.deepcopy(value_to_apply)
        else:
            base_dictionary[key_to_apply] = copy.deepcopy(value_to_apply)

    return base_dictionary


def _parse_inputconfig(input_config, yaml):
    """
    Utility function to parse input config into a dictionary.
    """
    if isinstance(input_config, dict):
        return input_config

    if isinstance(input_config, str):
        input_config = Path(Path.cwd(), input_config)

    with open(input_config) as fp:
        return yaml.load(fp)


def create_test_config(default_config, diff_config, output_path):
    """
    This function takes as input a default dictionary-like object and a
    difference dictionary-like object, merges the changes of the latter into
    the former, and saves the output in the desired output path.

    Inputs
    ------
    default_config : dict or path-like
        Default dictionary-like object provided as
        a dictionary or a path to a yaml file
    diff_config : dict or path-like
        Difference dictionary-like object provided as
        a dictionary or a path to a yaml file
    output_path : path-like
        Output path where the merged dictionary is saved

    Outputs
    -------
    - merged dictionary
    """

    # Load yaml files
    yaml = YAML()

    default_config = _parse_inputconfig(default_config, yaml)
    diff_config = _parse_inputconfig(diff_config, yaml)

    # create updated yaml
    merged_config = update(copy.deepcopy(default_config), diff_config)

    # Output path
    if isinstance(output_path, str):
        output_path = Path(Path.cwd(), output_path)

    # Save file
    yaml.dump(merged_config, output_path)

    return merged_config


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_test_configs")

    # Input paths
    fp_baseconfig = snakemake.input.base_config
    fp_update_file_list = snakemake.input.update_file_list
    fp_output_file_list = snakemake.output.tmp_test_configs

    for finput, foutput in zip(fp_update_file_list, fp_output_file_list):
        create_test_config(fp_baseconfig, finput, foutput)
