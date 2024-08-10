# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import sys

import pypsa

sys.path.append("./scripts")

from prepare_network import add_co2limit, download_emission_data, emission_extractor

emissions_file_url = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/v60_GHG/CO2_excl_short-cycle_org_C/v60_GHG_CO2_excl_short-cycle_org_C_1970_2018.zip"
emissions_file_name = "v60_CO2_excl_short-cycle_org_C_1970_2018.xls"
emissions_sheet_name = "v6.0_EM_CO2_fossil_IPCC1996"
automatic_emission_base_year = 1990
country_names = ["DE", "IT", "NG"]
co2limit = 7.75e7


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
    assert output_series.index.tolist() == ["NGA", "DEU", "ITA"]
    assert output_series.values.tolist() == [
        5698.76187,
        381475.887377666,
        123981.6458729,
    ]


def test_add_co2limit():
    test_network_de = pypsa.examples.scigrid_de(from_master=True)
    number_years = test_network_de.snapshot_weightings.objective.sum() / 8760.0
    add_co2limit(test_network_de, co2limit, number_years)
    assert (
        test_network_de.global_constraints.carrier_attribute.values[0]
        == "co2_emissions"
    )
    assert test_network_de.global_constraints.sense.values[0] == "<="
    assert test_network_de.global_constraints.constant.values[0] == 212328.76712328766
