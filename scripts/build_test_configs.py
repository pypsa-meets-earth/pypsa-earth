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
from pathlib import Path
from typing import Any

from ruamel.yaml import YAML


def update(
    d: dict[str, Any],
    u: collections.abc.Mapping[str, Any],
) -> dict[str, Any]:
    """
    Recursively merge mappings into a dictionary.

    Nested mappings in ``u`` are merged into the corresponding keys of ``d``;
    all other values in ``u`` overwrite values in ``d``.

    Parameters
    ----------
    d : dict
        Base dictionary to update.
    u : mapping
        Mapping whose items are merged into ``d``.

    Returns
    -------
    dict
        Updated dictionary (same object as ``d``).
    """
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def _parse_inputconfig(
    input_config: dict[str, Any] | str | Path,
    yaml: YAML,
) -> dict[str, Any]:
    """
    Parse a configuration object or YAML file into a dictionary.

    Parameters
    ----------
    input_config : dict or path-like
        Configuration provided as a dictionary or a path to a YAML file.
    yaml : ruamel.yaml.YAML
        YAML parser/loader instance.

    Returns
    -------
    dict
        Parsed configuration dictionary.
    """
    if isinstance(input_config, dict):
        return input_config

    if isinstance(input_config, str):
        input_config = Path(Path.cwd(), input_config)

    with open(input_config) as fp:
        return yaml.load(fp)


def create_test_config(
    default_config: dict[str, Any] | str | Path,
    diff_config: dict[str, Any] | str | Path,
    output_path: str | Path,
) -> dict[str, Any]:
    """
    Merge a default config with a diff config and write the result to disk.

    Takes a default dictionary-like object and a difference dictionary-like
    object, merges the changes of the latter into the former, and saves the
    output at the desired path.

    Parameters
    ----------
    default_config : dict or path-like
        Default configuration provided as a dictionary or a path to a YAML file.
    diff_config : dict or path-like
        Difference configuration provided as a dictionary or a path to a YAML
        file. Values in this config override or extend ``default_config``.
    output_path : path-like
        Output path where the merged configuration is saved.

    Returns
    -------
    dict
        Merged configuration dictionary.
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
