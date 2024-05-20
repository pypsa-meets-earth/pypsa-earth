# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import pathlib
import numpy as np
from scripts._helpers import (
    country_name_2_two_digits,
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
