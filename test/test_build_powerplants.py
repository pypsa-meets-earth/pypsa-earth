# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import pathlib
import sys

import pandas as pd
import yaml

sys.path.append("./scripts")

from test.conftest import get_config_dict

from build_powerplants import add_power_plants, replace_natural_gas_technology

path_cwd = pathlib.Path.cwd()


def test_replace_natural_gas_technology():
    """
    Verify what returned by replace_natural_gas_technology.
    """
    input_df = pd.DataFrame(
        {
            "Fueltype": [
                "Natural Gas",
                "Oil",
                "Natural Gas",
                "Natural Gas",
                "Natural Gas",
                "Natural Gas",
                "Natural Gas",
                "Natural Gas",
                "Natural Gas",
                "Natural Gas",
                "Natural Gas",
                "Natural Gas",
                "Hydro",
            ],
            "Technology": [
                "Steam Turbine",
                "Combustion Engine",
                "NG",
                "Ng",
                "NG/FO",
                "Ng/Fo",
                "NG/D",
                "LNG",
                "CCGT/D",
                "CCGT/FO",
                "LCCGT",
                "CCGT/Fo",
                "Reservoir",
            ],
        }
    )

    reference_df = pd.DataFrame(
        {
            "Fueltype": [
                "CCGT",
                "Oil",
                "CCGT",
                "CCGT",
                "OCGT",
                "OCGT",
                "OCGT",
                "OCGT",
                "CCGT",
                "CCGT",
                "CCGT",
                "CCGT",
                "Hydro",
            ],
            "Technology": [
                "CCGT",
                "Combustion Engine",
                "CCGT",
                "CCGT",
                "OCGT",
                "OCGT",
                "OCGT",
                "OCGT",
                "CCGT",
                "CCGT",
                "CCGT",
                "CCGT",
                "Reservoir",
            ],
        }
    )
    modified_df = replace_natural_gas_technology(input_df)
    comparison_df = modified_df.compare(reference_df)
    assert comparison_df.empty


def test_add_power_plants(get_config_dict):
    config_dict = get_config_dict
    custom_power_plants_file_path = pathlib.Path(
        path_cwd, "test", "test_data", "custom_NG_powerplants.csv"
    )
    pm_config_path = pathlib.Path(path_cwd, "configs", "powerplantmatching_config.yaml")
    with open(pm_config_path, "r") as f:
        power_plants_config = yaml.safe_load(f)
    ppl_query = config_dict["electricity"]["powerplants_filter"]

    config_dict["countries"] = ["NG"]

    # merge
    config_dict["electricity"]["custom_powerplants"] = "merge"
    custom_power_plant_query = config_dict["electricity"]["custom_powerplants"]
    if isinstance(ppl_query, str):
        power_plants_config["main_query"] = ppl_query
    countries_names = ["Nigeria"]
    power_plants_config["target_countries"] = countries_names
    ppl_merge = add_power_plants(
        custom_power_plants_file_path,
        power_plants_config,
        custom_power_plant_query,
        countries_names,
    )
    assert ppl_merge.shape == (64, 20)

    # replace
    config_dict["electricity"]["custom_powerplants"] = "replace"
    custom_power_plant_query = config_dict["electricity"]["custom_powerplants"]
    if isinstance(ppl_query, str):
        power_plants_config["main_query"] = ppl_query
    countries_names = ["Nigeria"]
    power_plants_config["target_countries"] = countries_names
    ppl_merge = add_power_plants(
        custom_power_plants_file_path,
        power_plants_config,
        custom_power_plant_query,
        countries_names,
    )
    assert ppl_merge.shape == (33, 19)

    # false
    config_dict["electricity"]["custom_powerplants"] = "false"
    custom_power_plant_query = config_dict["electricity"]["custom_powerplants"]
    if isinstance(ppl_query, str):
        power_plants_config["main_query"] = ppl_query
    countries_names = ["Nigeria"]
    power_plants_config["target_countries"] = countries_names
    ppl_merge = add_power_plants(
        custom_power_plants_file_path,
        power_plants_config,
        custom_power_plant_query,
        countries_names,
    )
    assert ppl_merge.shape == (31, 18)
