# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-


import sys

from pandas import Timestamp

sys.path.append("./scripts")

from test.conftest import get_power_network_ac_dc_meshed, get_power_network_scigrid_de

from prepare_network import (
    add_co2limit,
    add_emission_prices,
    add_gaslimit,
    average_every_nhours,
    download_emission_data,
    emission_extractor,
    enforce_autarky,
    set_line_nom_max,
    set_line_s_max_pu,
)

automatic_emission_base_year = 1990
country_names = ["DE", "IT", "NG"]
co2limit = 7.75e7


def test_download_emission_data():
    """
    Verify what returned by download_emission_data.
    """
    file_name = download_emission_data()
    assert file_name == "v60_CO2_excl_short-cycle_org_C_1970_2018.xls"


def test_emission_extractor():
    """
    Verify what returned by emission_extractor.
    """
    output_series = emission_extractor(
        "v60_CO2_excl_short-cycle_org_C_1970_2018.xls",
        automatic_emission_base_year,
        country_names,
    )
    assert output_series.index.tolist() == ["NGA", "DEU", "ITA"]
    assert output_series.values.tolist() == [
        5698.76187,
        381475.887377666,
        123981.6458729,
    ]


def test_add_co2limit(get_power_network_scigrid_de):
    """
    Verify what returned by add_co2limit.
    """
    test_network_de = get_power_network_scigrid_de
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


def test_add_gaslimit(get_power_network_scigrid_de):
    """
    Verify what returned by add_gaslimit.
    """
    test_network_de = get_power_network_scigrid_de
    test_network_de.add("Carrier", "OCGT")
    number_years = test_network_de.snapshot_weightings.objective.sum() / 8760.0
    add_gaslimit(test_network_de, number_years, number_years)
    assert test_network_de.global_constraints.carrier_attribute.values[0] == "gas_usage"
    assert test_network_de.global_constraints.sense.values[0] == "<="
    assert (
        test_network_de.global_constraints.constant.values[0]
        == number_years * number_years
    )


def test_add_emission_prices(get_power_network_ac_dc_meshed):
    """
    Verify what returned by add_emission_prices.
    """
    test_network_ac_dc_meshed = get_power_network_ac_dc_meshed
    add_emission_prices(
        test_network_ac_dc_meshed, emission_prices={"co2": 1.0}, exclude_co2=False
    )
    assert test_network_ac_dc_meshed.generators["marginal_cost"].index.tolist() == [
        "Manchester Wind",
        "Manchester Gas",
        "Norway Wind",
        "Norway Gas",
        "Frankfurt Wind",
        "Frankfurt Gas",
    ]
    assert test_network_ac_dc_meshed.generators["marginal_cost"].values.tolist() == [
        0.11,
        5.218030132047726,
        0.09,
        6.565421697244602,
        0.1,
        4.768788024122113,
    ]


def test_set_line_s_max_pu(get_power_network_scigrid_de):
    """
    Verify what returned by set_line_s_max_pu.
    """
    test_network_de = get_power_network_scigrid_de
    s_max_pu_new_value = 3.0
    set_line_s_max_pu(test_network_de, s_max_pu_new_value)
    assert test_network_de.lines["s_max_pu"].unique()[0] == s_max_pu_new_value


def test_average_every_nhours(get_power_network_scigrid_de):
    """
    Verify what returned by average_every_nhours.
    """
    test_network_de = get_power_network_scigrid_de

    # The input network is already sampled in 1H snapshots.
    # Hence, average_every_nhours should not change anything
    new_network_1h = average_every_nhours(test_network_de, "1H")
    assert test_network_de.snapshots.tolist() == new_network_1h.snapshots.tolist()

    # Re-sample to 4H intervals
    new_network_4h = average_every_nhours(test_network_de, "4H")
    assert new_network_4h.snapshots.tolist() == [
        Timestamp("2011-01-01 00:00:00"),
        Timestamp("2011-01-01 04:00:00"),
        Timestamp("2011-01-01 08:00:00"),
        Timestamp("2011-01-01 12:00:00"),
        Timestamp("2011-01-01 16:00:00"),
        Timestamp("2011-01-01 20:00:00"),
    ]


def test_enforce_autarky_only_crossborder_false(get_power_network_ac_dc_meshed):
    """
    Verify what returned by enforce_autarky when only_crossborder is False.
    """
    # --> it removes all lines and all DC links
    test_network_no_cross_border = get_power_network_ac_dc_meshed

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


def test_enforce_autarky_only_crossborder_true(get_power_network_ac_dc_meshed):
    """
    Verify what returned by enforce_autarky when only_crossborder is True.
    """
    # --> it removes links and lines that cross borders
    test_network_with_cross_border = get_power_network_ac_dc_meshed
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
    print(output_component_dict_with_cross_border)

    assert (
        output_component_dict_with_cross_border
        == reference_component_dict_with_cross_border
    )


def test_set_line_nom_max(get_power_network_ac_dc_meshed):
    """
    Verify what returned by set_line_nom_max.
    """
    test_network = get_power_network_ac_dc_meshed
    s_nom_max_value = 5.0
    p_nom_max_value = 10.0
    set_line_nom_max(
        test_network, s_nom_max_set=s_nom_max_value, p_nom_max_set=p_nom_max_value
    )
    assert test_network.lines.s_nom_max.unique()[0] == s_nom_max_value
    assert test_network.links.p_nom_max.unique()[0] == p_nom_max_value
