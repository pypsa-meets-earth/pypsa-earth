# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import os
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.append("./scripts")

from build_demand_profiles import get_gegis_regions, get_load_paths_gegis

path_cwd = pathlib.Path.cwd()


def test_get_gegis_regions():
    """
    Verify what returned by get_gegis_regions.
    """
    output_regions = get_gegis_regions(["NG", "IT"])
    assert output_regions == ["Africa", "Europe"]


def test_get_load_paths_gegis():
    """
    Verify what returned by get_load_paths_gegis.
    """
    config = {
        "countries": ["NG", "IT"],
        "load_options": {
            "ssp": "ssp2-2.6",
            "weather_year": 2013,
            "prediction_year": 2030,
            "scale": 1,
        },
    }
    load_data_paths = get_load_paths_gegis("data", config)
    reference_list = [
        os.path.join("data", "ssp2-2.6", "2030", "era5_2013", "Africa.nc"),
        os.path.join("data", "ssp2-2.6", "2030", "era5_2013", "Europe.nc"),
    ]
    assert load_data_paths == reference_list
