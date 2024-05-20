# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import pathlib
import numpy as np
from scripts._helpers import (
    country_name_2_two_digits,
    get_conv_factors,
    safe_divide,
    sets_path_to_root,
    two_2_three_digits_country,
    two_digits_2_name_country,
    three_2_two_digits_country,)

path_cwd = pathlib.Path.cwd()

def test_sets_path_to_root():
    """
    Verify that the root folder returned by sets_path_to_root is pypsa-earth
    """
    sets_path_to_root("pypsa-earth")
    current_path = pathlib.Path.cwd()
    assert current_path == path_cwd


def test_two_2_three_digits_country():
    """
    Verify the conversion from two-digit to three-digit country code
    """
    # Afghanistan
    assert two_2_three_digits_country("AF") == "AFG"
    # American Samoa
    assert two_2_three_digits_country("AS") == "ASM"
    # Aruba
    assert two_2_three_digits_country("AW") == "ABW"
    # Germany
    assert two_2_three_digits_country("DE") == "DEU"
    # Micronesia (Federated States of)
    assert two_2_three_digits_country("FM") == "FSM"


def test_three_2_two_digits_country():
    """
    Verify the conversion from three-digit to two-digit country code
    """
    # Afghanistan
    assert "AF" == three_2_two_digits_country("AFG")
    # American Samoa
    assert "AS" == three_2_two_digits_country("ASM")
    # Aruba
    assert "AW" == three_2_two_digits_country("ABW")
    # Germany
    assert "DE" == three_2_two_digits_country("DEU")
    # Micronesia (Federated States of)
    assert "FM" == three_2_two_digits_country("FSM")


def test_two_digits_2_name_country():
    """
    Verify the conversion from two-digit country code to country name
    """
    # Micronesia (Federated States of)
    assert "Micronesia, Fed. Sts." == two_digits_2_name_country("FM")
    assert "Fed. Sts. Micronesia" == two_digits_2_name_country("FM", nocomma=True)
    assert "Sts. Micronesia" == two_digits_2_name_country("FM", nocomma=True, remove_start_words=["Fed. "])


def test_country_name_2_two_digits():
    """
    Verify the conversion from country name to two-digit country code
    """
    # Afghanistan
    assert "AF" == country_name_2_two_digits("Afghanistan")
    # American Samoa
    assert "AS" == country_name_2_two_digits("American Samoa")
    # Aruba
    assert "AW" == country_name_2_two_digits("Aruba")
    # Germany
    assert "DE" == country_name_2_two_digits("Germany")
    # Micronesia (Federated States of)
    assert "FM" == country_name_2_two_digits("Micronesia")


def test_safe_divide():
    """
    Verify that the method safe_divide prevents divisions by vanishing denominator
    """
    assert safe_divide(3.0,2.0) == 1.5
    assert np.isnan(safe_divide(3.0, 0.0))


def test_get_conv_factors():
    """
    Verify that the conversion factors returned by get_conv_factors are correct
    """
    conversion_factors_dict = get_conv_factors("industry")
    assert conversion_factors_dict["additives and oxygenates"] == 0.008333
    assert conversion_factors_dict["anthracite"] == 0.005
    assert conversion_factors_dict["aviation gasoline"] == 0.01230
    assert conversion_factors_dict["bagasse"] == 0.002144
    assert conversion_factors_dict["biodiesel"] == 0.01022
    assert conversion_factors_dict["biogasoline"] == 0.007444
    assert conversion_factors_dict["bio jet kerosene"] == 0.011111
    assert conversion_factors_dict["bitumen"] == 0.01117
    assert conversion_factors_dict["brown coal"] == 0.003889
    assert conversion_factors_dict["brown coal briquettes"] == 0.00575
    assert conversion_factors_dict["charcoal"] == 0.00819
    assert conversion_factors_dict["coal tar"] == 0.007778
    assert conversion_factors_dict["coke oven coke"] == 0.0078334
    assert conversion_factors_dict["coking coal"] == 0.007833
    assert conversion_factors_dict["conventional crude oil"] == 0.01175
    assert conversion_factors_dict["crude petroleum"] == 0.011750
    assert conversion_factors_dict["fuel oil"] == 0.01122
    assert conversion_factors_dict["fuelwood"] == 0.00254
    assert conversion_factors_dict["gas coke"] == 0.007326
    assert conversion_factors_dict["gas oil / diesel oil"] == 0.01194
    assert conversion_factors_dict["gasoline type jet fuel"] == 0.01230
    assert conversion_factors_dict["hard coal"] == 0.007167
    assert conversion_factors_dict["kerosene type jet fuel"] == 0.01225
    assert conversion_factors_dict["lignite"] == 0.003889
    assert conversion_factors_dict["liquefied petroleum gas (LPG)"] == 0.01313
    assert conversion_factors_dict["lubricants"] == 0.01117
    assert conversion_factors_dict["motor gasoline"] == 0.01230
    assert conversion_factors_dict["naphtha"] == 0.01236
    assert conversion_factors_dict["natural gas liquids"] == 0.01228
    assert conversion_factors_dict["other bituminous coal"] == 0.005556
    assert conversion_factors_dict["patent fuel"] == 0.00575
    assert conversion_factors_dict["peat"] == 0.00271
    assert conversion_factors_dict["peat products"] == 0.00271
    assert conversion_factors_dict["petroleum coke"] == 0.009028
    assert conversion_factors_dict["refinery gas"] == 0.01375
    assert conversion_factors_dict["sub bituminous coal"] == 0.005555
    assert np.isnan(get_conv_factors("non-industry"))
