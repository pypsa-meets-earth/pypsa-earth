# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import sys

sys.path.append("./scripts")

from prepare_network import download_emission_data


def test_download_emission_data():
    emissions_file_url = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/v60_GHG/CO2_excl_short-cycle_org_C/v60_GHG_CO2_excl_short-cycle_org_C_1970_2018.zip"
    emissions_file_name = "v60_CO2_excl_short-cycle_org_C_1970_2018.xls"
    filename = download_emission_data(emissions_file_url, emissions_file_name)
    assert filename == emissions_file_name
