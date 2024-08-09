# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import sys

import numpy as np
import pandas as pd

sys.path.append("./scripts")

from prepare_network import download_emission_data, emission_extractor

emissions_file_url = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/v60_GHG/CO2_excl_short-cycle_org_C/v60_GHG_CO2_excl_short-cycle_org_C_1970_2018.zip"
emissions_file_name = "v60_CO2_excl_short-cycle_org_C_1970_2018.xls"
emissions_sheet_name = "v6.0_EM_CO2_fossil_IPCC1996"
automatic_emission_base_year = 1990
country_names = ["DE", "IT", "NG"]


def test_download_emission_data():
    filename = download_emission_data(emissions_file_url, emissions_file_name)
    assert filename == emissions_file_name


def test_emission_extractor():
    output_series = emission_extractor(
        emissions_file_name,
        emissions_sheet_name,
        automatic_emission_base_year,
        country_names,
    )
    reference_df = pd.DataFrame(
        {
            "Country_code_A3": ["NGA", "DEU", "ITA"],
            "emissions": [5698.761870, 381475.887378, 123981.645873],
        }
    )
    output_df = pd.DataFrame(
        {"Country_code_A3": output_series.index, "emissions": output_series.values}
    )
    print(reference_df)
    print(output_df)
    comparison_df = output_df.compare(reference_df)
    print(comparison_df)
    print(type(output_series.values))
    print(type(output_series.values[0]))
    print(type(output_series.values[1]))
    print(type(output_series.values[2]))
    print(type(5698.761870))
    assert comparison_df.empty
    # assert False
