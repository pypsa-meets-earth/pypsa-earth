# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Checks relevance of `config.default.yaml`.

Relevant Settings
-----------------

.. code:: yaml

.. seealso::
    Documentation of management of the configuration file

Inputs
------
- ``configs``
- ``config.default.yaml``

Description
-----------
Checks that all parameters of `config.default.yaml` are included into
configuration files in `config` folder
"""
import glob
import os
import sys

import yaml

sys.path.append("./scripts")

from _helpers import (
    compare_configs,
    configure_logging,
    create_logger,
)

logger = create_logger(__name__)


def read_yaml_file(filename):
    with open(filename, "r") as stream:
        res = yaml.safe_load(stream)
    return res


def load_yaml_files(directory, pattern="*.yaml"):
    merged_configs = {}

    for filepath in glob.glob(os.path.join(directory, pattern)):
        with open(filepath, "r") as file:
            current_config = yaml.safe_load(file)

            # dictionary keys should be unique in Python
            # while it can happen that entries of the same key
            # will be distributed across dirrerent configs/*.config.yaml
            duplicated_keys = current_config.keys() & merged_configs.keys()
            if duplicated_keys:
                logger.error(
                    f"There is duplication of a key across the config files for {duplicated_keys}"
                    + " That leads to overwriting and a potential loss of data"
                )
            merged_configs.update(current_config)

    return merged_configs


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("check_default_config")

    configure_logging(snakemake)

    configs_folder = snakemake.input["config_folder"]
    configs_yaml = load_yaml_files(directory=configs_folder, pattern="*.yaml")

    for file in ["config.default.yaml"]:
        default_config_yaml = read_yaml_file(file)

    config_compar_dict = compare_configs(
        current_config=default_config_yaml,
        benchmark_config=configs_yaml,
    )
