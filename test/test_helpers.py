# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import os
import pathlib
import shutil
import sys
from test.conftest import (
    _content_temp_file,
    _name_temp_file,
    _sub_temp_content_dir,
    _temp_content_dir,
    get_temp_file,
)

import fiona
import numpy as np
import pandas as pd

sys.path.append("./scripts")

from _helpers import (
    aggregate_fuels,
    country_name_2_two_digits,
    get_conv_factors,
    safe_divide,
    three_2_two_digits_country,
    two_2_three_digits_country,
    two_digits_2_name_country,
)

path_cwd = str(pathlib.Path.cwd())


original_commodity_data = [
    "Biogases",
    "Fuelwood",
    "of which: fishing",
    "Natural gas liquids",
    "Naphtha",
    "Motor Gasoline",
    "Motor gasoline",
    "Gasoline-type jet fuel",
    "Peat products",
    "Peat Products",
    "Direct use of geothermal heat",
    "Additives and Oxygenates",
    "Electricity",
    "Animal waste",
    "animal waste",
    "Refinery gas",
    "Refinery Gas",
    "Fuel oil",
    "Oil shale",
    "Oil Shale",
    "Lignite",
    "Falling water",
    "Petroleum coke",
    "Petroleum Coke",
    "Aviation gasoline",
    "Ethane",
    "Natural gas (including LNG)",
    "Natural gas",
    "Natural Gas (including LNG)",
    "Other bituminous coal",
    "Paraffin waxes",
    "Hard coal",
    "Coal",
    "Hrad coal",
    "Coke Oven Gas",
    "Gasworks Gas",
    "Brown coal briquettes",
    "Brown Coal Briquettes",
    "Liquefied petroleum gas (LPG)",
    "Liquified Petroleum Gas (LPG)",
    "Sub-bituminous coal",
    "Kerosene-type Jet Fuel",
    "Charcoal",
    "Heat",
    "Gas coke",
    "Gas Coke",
    "Patent fuel",
    "Peat (for fuel use)",
    "Peat",
    "Coal Tar",
    "Biogasoline",
    "Coking coal",
    "Electricity generating capacity",
    "Anthracite",
    "Coke oven coke",
    "Coke-oven coke",
    "Coke Oven Coke",
    "Conventional crude oil",
    "Crude petroleum",
    "Brown coal",
    "Lignite brown coal",
    "Lignite brown coal- recoverable resources",
    "Biodiesel",
    "Lubricants",
    "Black Liquor",
    "Gas Oil/ Diesel Oil",
    "Gas Oil/ Diesel Oil ",
    "Gas Oil/Diesel Oil",
    "Bagasse",
    "Direct use of solar thermal heat",
    "Bio jet kerosene",
    "Blast Furnace Gas",
    "Blast furnace gas",
    "Bitumen",
]

modified_commodity_data = [
    "biogases",
    "fuelwood",
    "of which: fishing",
    "natural gas liquids",
    "naphtha",
    "motor gasoline",
    "gasoline-type jet fuel",
    "peat products",
    "direct use of geothermal heat",
    "additives and oxygenates",
    "electricity",
    "animal waste",
    "refinery gas",
    "fuel oil",
    "oil shale",
    "lignite",
    "falling water",
    "petroleum coke",
    "aviation gasoline",
    "ethane",
    "natural gas (including lng)",
    "natural gas",
    "other bituminous coal",
    "paraffin waxes",
    "hard coal",
    "coal",
    "coke-oven gas",
    "gasworks gas",
    "brown coal briquettes",
    "liquefied petroleum gas (lpg)",
    "sub-bituminous coal",
    "kerosene-type jet fuel",
    "charcoal",
    "heat",
    "gas coke",
    "patent fuel",
    "peat (for fuel use)",
    "peat",
    "coal tar",
    "biogasoline",
    "coking coal",
    "electricity generating capacity",
    "anthracite",
    "coke-oven coke",
    "conventional crude oil",
    "crude petroleum",
    "brown coal",
    "lignite brown coal",
    "lignite brown coal - recoverable resources",
    "biodiesel",
    "lubricants",
    "black liquor",
    "gas oil/ diesel oil",
    "bagasse",
    "direct use of solar thermal heat",
    "bio jet kerosene",
    "blast furnace gas",
    "bitumen",
]

original_commodity_dataframe = pd.DataFrame(
    original_commodity_data, columns=["Commodity"]
)
modified_commodity_dataframe = pd.DataFrame(
    modified_commodity_data, columns=["Commodity"]
)


# def test_build_directory(get_temp_folder, tmpdir):
#     """
#     Verify the directory tree returned by build_directory()

#     Please note:
#     -) build_directory(path, just_parent_directory=True) is equivalent to os.makedirs(os.path.dirname(path)).
#     Given a path tmpdir/temp_content_dir/sub_temp_content_dir, it will create just tmpdir/temp_content_dir/
#     -) build_directory(path, just_parent_directory=False) is equivalent to os.makedirs(path). Given a path
#     tmpdir/temp_content_dir/sub_temp_content_dir, it will create tmpdir/temp_content_dir/sub_temp_content_dir
#     """

#     # test with pathlib
#     build_directory(get_temp_folder, just_parent_directory=True)
#     just_parent_list_pathlib = []
#     for root, dirs, files in os.walk(tmpdir):
#         just_parent_list_pathlib.append(str(os.path.join(root)))

#     assert len(just_parent_list_pathlib) == 2
#     assert just_parent_list_pathlib[0] == str(tmpdir)
#     assert just_parent_list_pathlib[1] == str(tmpdir.join(_temp_content_dir))

#     # remove the temporary folder tmpdir/temp_content_dir/
#     shutil.rmtree(pathlib.Path(tmpdir, _temp_content_dir))

#     # test with os.makedirs. Please note for exist_ok=False,
#     # a FileExistsError is raised if the target directory
#     # already exists. Hence, setting exist_ok=False ensures
#     # that the removal with shutil.rmtree was successful
#     os.makedirs(os.path.dirname(get_temp_folder), exist_ok=False)
#     just_parent_list_os = []
#     for root, dirs, files in os.walk(tmpdir):
#         just_parent_list_os.append(str(os.path.join(root)))

#     assert just_parent_list_pathlib == just_parent_list_os

#     # test with pathlib
#     build_directory(get_temp_folder, just_parent_directory=False)
#     full_tree_list_pathlib = []
#     for root, dirs, files in os.walk(tmpdir):
#         full_tree_list_pathlib.append(str(os.path.join(root)))

#     assert len(full_tree_list_pathlib) == 3
#     assert full_tree_list_pathlib[0] == str(tmpdir)
#     assert full_tree_list_pathlib[1] == str(tmpdir.join(_temp_content_dir))
#     assert full_tree_list_pathlib[2] == str(
#         tmpdir.join(_temp_content_dir, _sub_temp_content_dir)
#     )

#     # remove the temporary folder tmpdir/temp_content_dir/*
#     shutil.rmtree(pathlib.Path(tmpdir, _temp_content_dir))

#     # test with os.makedirs. Please note for exist_ok=False,
#     # a FileExistsError is raised if the target directory
#     # already exists. Hence, setting exist_ok=False ensures
#     # that the removal with shutil.rmtree was successful
#     os.makedirs(get_temp_folder, exist_ok=False)
#     full_tree_list_os = []
#     for root, dirs, files in os.walk(tmpdir):
#         full_tree_list_os.append(str(os.path.join(root)))

#     assert full_tree_list_os == full_tree_list_pathlib


def test_two_2_three_digits_country():
    """
    Verify the conversion from two-digit to three-digit country code.
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
    Verify the conversion from three-digit to two-digit country code.
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
    Verify the conversion from two-digit country code to country name.
    """
    # Micronesia (Federated States of)
    assert "Micronesia, Fed. Sts." == two_digits_2_name_country("FM")
    assert "Fed. Sts." == two_digits_2_name_country(
        "FM", remove_start_words=["Micronesia, "]
    )
    # Democratic Republic of the Congo
    assert "DR Congo" == two_digits_2_name_country("CD")


def test_country_name_2_two_digits():
    """
    Verify the conversion from country name to two-digit country code.
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
    Verify that the method safe_divide prevents divisions by vanishing
    denominator.
    """
    numerator = pd.Series([1.0, 2.0], index=[0, 1])
    denominator = 2.0
    assert (safe_divide(numerator, denominator) == pd.Series([0.5, 1.0])).all()


def test_get_conv_factors():
    """
    Verify that the conversion factors returned by get_conv_factors are
    correct.
    """
    conversion_factors_dict = get_conv_factors("industry")
    assert conversion_factors_dict["Additives and Oxygenates"] == 0.008333
    assert conversion_factors_dict["Oil shale"] == 0.00247
