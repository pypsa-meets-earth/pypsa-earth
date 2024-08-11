# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import sys

import pypsa

sys.path.append("./scripts")

from prepare_network import (
    add_co2limit,
    add_emission_prices,
    add_gaslimit,
    download_emission_data,
    emission_extractor,
    enforce_autarky,
    set_line_s_max_pu,
)

emissions_file_url = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/v60_GHG/CO2_excl_short-cycle_org_C/v60_GHG_CO2_excl_short-cycle_org_C_1970_2018.zip"
emissions_file_name = "v60_CO2_excl_short-cycle_org_C_1970_2018.xls"
emissions_sheet_name = "v6.0_EM_CO2_fossil_IPCC1996"
automatic_emission_base_year = 1990
country_names = ["DE", "IT", "NG"]
co2limit = 7.75e7


def test_download_emission_data():
    file_name = download_emission_data(emissions_file_url, emissions_file_name)
    assert file_name == emissions_file_name


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
    assert (
        test_network_de.global_constraints.constant.values[0] == co2limit * number_years
    )


def test_add_gaslimit():
    test_network_de = pypsa.examples.scigrid_de(from_master=True)
    test_network_de.add("Carrier", "OCGT")
    number_years = test_network_de.snapshot_weightings.objective.sum() / 8760.0
    add_gaslimit(test_network_de, number_years, number_years)
    assert test_network_de.global_constraints.carrier_attribute.values[0] == "gas_usage"
    assert test_network_de.global_constraints.sense.values[0] == "<="
    assert (
        test_network_de.global_constraints.constant.values[0]
        == number_years * number_years
    )


def test_set_line_s_max_pu():
    test_network_de = pypsa.examples.scigrid_de(from_master=True)
    s_max_pu_new_value = 3.0
    set_line_s_max_pu(test_network_de, s_max_pu_new_value)
    assert test_network_de.lines["s_max_pu"].unique()[0] == s_max_pu_new_value


def test_enforce_autarky():

    # test with only_crossborder=False
    test_network_no_cross_border = pypsa.examples.ac_dc_meshed(from_master=True)

    bus_country_list = ["UK", "UK", "UK", "UK", "DE", "DE", "DE", "NO", "NO"]
    test_network_no_cross_border.buses["country"] = bus_country_list
    test_network_no_cross_border.links["carrier"] = "DC"

    enforce_autarky(test_network_no_cross_border, only_crossborder=False)

    output_component_dict_no_cross_border = {}
    for c in test_network_no_cross_border.iterate_components(
        list(test_network_no_cross_border.components.keys())[2:]
    ):
        output_component_dict_no_cross_border[c.name] = len(c.df)

    reference_component_dict_no_cross_border = {
        "Bus": 9,
        "Carrier": 3,
        "GlobalConstraint": 1,
        "LineType": 34,
        "TransformerType": 14,
        "Load": 6,
        "Generator": 6,
    }
    assert (
        output_component_dict_no_cross_border
        == reference_component_dict_no_cross_border
    )

    # test with only_crossborder=True
    test_network_with_cross_border = pypsa.examples.ac_dc_meshed(from_master=True)
    bus_country_list = ["UK", "UK", "UK", "UK", "DE", "DE", "DE", "NO", "NO"]
    test_network_with_cross_border.buses["country"] = bus_country_list
    test_network_with_cross_border.links["carrier"] = "DC"

    enforce_autarky(test_network_with_cross_border, only_crossborder=True)

    output_component_dict_with_cross_border = {}
    for c in test_network_with_cross_border.iterate_components(
        list(test_network_with_cross_border.components.keys())[2:]
    ):
        output_component_dict_with_cross_border[c.name] = len(c.df)

    reference_component_dict_with_cross_border = {
        "Bus": 9,
        "Carrier": 3,
        "GlobalConstraint": 1,
        "Line": 4,
        "LineType": 34,
        "TransformerType": 14,
        "Link": 3,
        "Load": 6,
        "Generator": 6,
    }
    assert (
        output_component_dict_with_cross_border
        == reference_component_dict_with_cross_border
    )


# I do not understand this method!
# def test_add_emission_prices():
#     test_network_de = pypsa.examples.scigrid_de(from_master=True)
#     test_network_de.add("Carrier", "OCGT_emissions", co2_emissions=2.0)
#     test_network_de.add("Carrier", "coal_emissions", co2_emissions=2.0)
#     test_network_de.add("Generator", "OCGT_emissions")
#     test_network_de.add("Generator", "coal_emissions")
#     ep = (
#         pd.Series({"co2": 1.0}).rename(lambda x: x + "_emissions")
#         * test_network_de.carriers.filter(like="_emissions")
#     ).sum(axis=1)
#     print("first term", test_network_de.generators.carrier.map(ep))
#     print("second term", test_network_de.generators.efficiency)
#     #print("before", test_network_de.generators["marginal_cost"])
#     #print("=====")
#     add_emission_prices(test_network_de, emission_prices={"co2": 1.0}, exclude_co2=False)
#     #print("after", test_network_de.generators["marginal_cost"])
#     assert False
