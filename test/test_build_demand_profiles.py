# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import pathlib
import sys

sys.path.append("./scripts")

from test.conftest import get_config_dict

from _helpers import create_country_list, get_path
from build_demand_profiles import get_gegis_regions, get_load_paths_gegis

path_cwd = pathlib.Path.cwd()


def test_get_gegis_regions():
    """
    Verify what returned by get_gegis_regions.
    """
    output_regions = get_gegis_regions(["NG", "IT"])
    assert output_regions == ["Africa", "Europe"]


def test_get_load_paths_gegis(get_config_dict):
    """
    Verify what returned by get_load_paths_gegis.
    """
    config_dict = get_config_dict
    load_data_paths = get_load_paths_gegis("data", config_dict)
    reference_list = [
        get_path("data", "ssp2-2.6", "2030", "era5_2013", "Africa.nc"),
    ]
    print(load_data_paths)
    assert load_data_paths == reference_list
