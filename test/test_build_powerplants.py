# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import pathlib
import sys

import pandas as pd

sys.path.append("./scripts")

from test.conftest import get_config_dict

from _helpers import create_country_list, get_path
from build_powerplants import replace_natural_gas_technology

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
