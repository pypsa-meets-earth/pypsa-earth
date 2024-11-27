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
    build_directory,
    change_to_script_dir,
    country_name_2_two_digits,
    download_gadm,
    get_base_dir,
    get_config_default_path,
    get_conv_factors,
    get_current_directory_path,
    get_gadm_filename,
    get_gadm_url,
    get_path,
    get_path_size,
    get_relative_path,
    modify_commodity,
    normed,
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


def test_build_directory(get_temp_folder, tmpdir):
    """
    Verify the directory tree returned by build_directory()

    Please note:
    -) build_directory(path, just_parent_directory=True) is equivalent to os.makedirs(os.path.dirname(path)).
    Given a path tmpdir/temp_content_dir/sub_temp_content_dir, it will create just tmpdir/temp_content_dir/
    -) build_directory(path, just_parent_directory=False) is equivalent to os.makedirs(path). Given a path
    tmpdir/temp_content_dir/sub_temp_content_dir, it will create tmpdir/temp_content_dir/sub_temp_content_dir
    """

    # test with pathlib
    build_directory(get_temp_folder, just_parent_directory=True)
    just_parent_list_pathlib = []
    for root, dirs, files in os.walk(tmpdir):
        just_parent_list_pathlib.append(str(get_path(root)))

    assert len(just_parent_list_pathlib) == 2
    assert just_parent_list_pathlib[0] == str(tmpdir)
    assert just_parent_list_pathlib[1] == str(tmpdir.join(_temp_content_dir))

    # remove the temporary folder tmpdir/temp_content_dir/
    shutil.rmtree(pathlib.Path(tmpdir, _temp_content_dir))

    # test with os.makedirs. Please note for exist_ok=False,
    # a FileExistsError is raised if the target directory
    # already exists. Hence, setting exist_ok=False ensures
    # that the removal with shutil.rmtree was successful
    os.makedirs(os.path.dirname(get_temp_folder), exist_ok=False)
    just_parent_list_os = []
    for root, dirs, files in os.walk(tmpdir):
        just_parent_list_os.append(str(get_path(root)))

    assert just_parent_list_pathlib == just_parent_list_os

    # test with pathlib
    build_directory(get_temp_folder, just_parent_directory=False)
    full_tree_list_pathlib = []
    for root, dirs, files in os.walk(tmpdir):
        full_tree_list_pathlib.append(str(get_path(root)))

    assert len(full_tree_list_pathlib) == 3
    assert full_tree_list_pathlib[0] == str(tmpdir)
    assert full_tree_list_pathlib[1] == str(tmpdir.join(_temp_content_dir))
    assert full_tree_list_pathlib[2] == str(
        tmpdir.join(_temp_content_dir, _sub_temp_content_dir)
    )

    # remove the temporary folder tmpdir/temp_content_dir/*
    shutil.rmtree(pathlib.Path(tmpdir, _temp_content_dir))

    # test with os.makedirs. Please note for exist_ok=False,
    # a FileExistsError is raised if the target directory
    # already exists. Hence, setting exist_ok=False ensures
    # that the removal with shutil.rmtree was successful
    os.makedirs(get_temp_folder, exist_ok=False)
    full_tree_list_os = []
    for root, dirs, files in os.walk(tmpdir):
        full_tree_list_os.append(str(get_path(root)))

    assert full_tree_list_os == full_tree_list_pathlib


def test_change_to_script_dir():
    """
    Verify the path returned by change_to_script_dir()
    """
    change_to_script_dir(__file__)
    assert str(pathlib.Path.cwd()) == path_cwd + os.sep + "test"
    change_to_script_dir(".")
    assert str(pathlib.Path.cwd()) == path_cwd


def test_get_path():
    """
    Verify the path returned by get_path()
    """
    file_name_path_one = get_path(
        path_cwd,
        "sub_path_1",
        "sub_path_2",
        "sub_path_3",
        "sub_path_4",
        "sub_path_5",
        "file.nc",
    )
    file_name_path_one_list_unpacked = get_path(
        path_cwd,
        *[
            "sub_path_1",
            "sub_path_2",
            "sub_path_3",
            "sub_path_4",
            "sub_path_5",
            "file.nc",
        ],
    )
    path_name_path_two = get_path(
        pathlib.Path(__file__).parent, "..", "logs", "rule.log"
    )
    assert str(file_name_path_one) == os.path.join(
        path_cwd,
        "sub_path_1",
        "sub_path_2",
        "sub_path_3",
        "sub_path_4",
        "sub_path_5",
        "file.nc",
    )
    assert str(file_name_path_one_list_unpacked) == os.path.join(
        path_cwd,
        "sub_path_1",
        "sub_path_2",
        "sub_path_3",
        "sub_path_4",
        "sub_path_5",
        "file.nc",
    )
    assert str(path_name_path_two) == str(
        pathlib.Path(__file__).parent.joinpath("..", "logs", "rule.log")
    )


def test_get_path_size(get_temp_file):
    """
    Verify the path size (in bytes) returned by get_path_size()
    """
    path = get_temp_file
    file_size = get_path_size(path)
    assert file_size == os.stat(path).st_size
    assert file_size == len(_content_temp_file)


def test_get_base_dir(get_temp_file):
    """
    Verify the base directory path returned by get_base_dir()
    """
    assert get_base_dir(get_temp_file) == os.path.abspath(os.path.join(os.path.dirname(get_temp_file), os.pardir))
    assert get_base_dir(__file__) == os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))


def test_get_config_default_path(get_temp_file):
    """
    Verify the base directory path returned by get_config_default_path()
    """
    base_dir = get_config_default_path(get_temp_file)
    assert get_config_default_path(base_dir) == str(get_path(base_dir, "config.default.yaml"))
    base_dir = get_config_default_path(__file__)
    assert get_config_default_path(base_dir) == str(get_path(base_dir, "config.default.yaml"))


def test_get_current_directory_path():
    """
    Verify the current directory path returned by get_current_directory_path()
    """
    path = get_current_directory_path()
    assert str(path) == os.getcwd()


def test_get_relative_path(get_temp_file):
    """
    Verify the relative path returned by get_relative_path()
    """
    path = get_temp_file
    # path relative to the parent directory of the temp file
    relative_path = get_relative_path(path, get_path(path).parent)
    assert str(relative_path) == _name_temp_file
    assert str(relative_path) == os.path.relpath(path, start=get_path(path).parent)


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
    assert "Federated States of Micronesia" == two_digits_2_name_country(
        "FM", name_string="name_official"
    )
    assert "States of Micronesia" == two_digits_2_name_country(
        "FM", name_string="name_official", remove_start_words=["Federated "]
    )
    # Democratic Republic of the Congo
    assert "DR Congo" == two_digits_2_name_country("CD")
    assert "Democratic Republic of the Congo" == two_digits_2_name_country(
        "CD", name_string="name_official"
    )
    assert "Republic of the Congo" == two_digits_2_name_country(
        "CD", name_string="name_official", remove_start_words=["Democratic "]
    )


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
    assert safe_divide(3.0, 2.0) == 1.5
    assert np.isnan(safe_divide(3.0, 0.0))


def test_get_conv_factors():
    """
    Verify that the conversion factors returned by get_conv_factors are
    correct.
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
    assert conversion_factors_dict["coke-oven coke"] == 0.0078334
    assert conversion_factors_dict["coke-oven gas"] == 0.000277
    assert conversion_factors_dict["coking coal"] == 0.007833
    assert conversion_factors_dict["conventional crude oil"] == 0.01175
    assert conversion_factors_dict["crude petroleum"] == 0.011750
    assert conversion_factors_dict["ethane"] == 0.01289
    assert conversion_factors_dict["fuel oil"] == 0.01122
    assert conversion_factors_dict["fuelwood"] == 0.00254
    assert conversion_factors_dict["gas coke"] == 0.007326
    assert conversion_factors_dict["gas oil/ diesel oil"] == 0.01194
    assert conversion_factors_dict["gasoline-type jet fuel"] == 0.01230
    assert conversion_factors_dict["hard coal"] == 0.007167
    assert conversion_factors_dict["kerosene-type jet fuel"] == 0.01225
    assert conversion_factors_dict["lignite"] == 0.003889
    assert conversion_factors_dict["liquefied petroleum gas (lpg)"] == 0.01313
    assert conversion_factors_dict["lubricants"] == 0.011166
    assert conversion_factors_dict["motor gasoline"] == 0.01230
    assert conversion_factors_dict["naphtha"] == 0.01236
    assert conversion_factors_dict["natural gas liquids"] == 0.01228
    assert conversion_factors_dict["oil shale"] == 0.00247
    assert conversion_factors_dict["other bituminous coal"] == 0.005556
    assert conversion_factors_dict["paraffin waxes"] == 0.01117
    assert conversion_factors_dict["patent fuel"] == 0.00575
    assert conversion_factors_dict["peat"] == 0.00271
    assert conversion_factors_dict["peat products"] == 0.00271
    assert conversion_factors_dict["petroleum coke"] == 0.009028
    assert conversion_factors_dict["refinery gas"] == 0.01375
    assert conversion_factors_dict["sub-bituminous coal"] == 0.005555
    assert np.isnan(get_conv_factors("non-industry"))


def test_modify_commodity():
    """
    Verify that modify_commodity returns the commodities in wished format.
    """
    new_commodity_dataframe = pd.DataFrame()
    new_commodity_dataframe["Commodity"] = (
        original_commodity_dataframe["Commodity"].map(modify_commodity).unique()
    )
    df = new_commodity_dataframe.compare(modified_commodity_dataframe)
    boolean_flag = df.empty
    if not boolean_flag:
        assert False


def test_aggregate_fuels():
    """
    Verify what is returned by aggregate_fuels.
    """
    assert np.isnan(aggregate_fuels("non-industry"))


def test_get_gadm_filename():
    """
    Verify what is returned by get_gadm_filename.
    """
    # Kosovo
    assert get_gadm_filename("XK") == "gadm41_XKO"
    # Clipperton island
    assert get_gadm_filename("CP") == "gadm41_XCL"
    # Saint-Martin
    assert get_gadm_filename("SX") == "gadm41_MAF"
    # French Southern Territories
    assert get_gadm_filename("TF") == "gadm41_ATF"
    # Aland
    assert get_gadm_filename("AX") == "gadm41_ALA"
    # British Indian Ocean Territory
    assert get_gadm_filename("IO") == "gadm41_IOT"
    # Cocos Islands
    assert get_gadm_filename("CC") == "gadm41_CCK"
    # Norfolk
    assert get_gadm_filename("NF") == "gadm41_NFK"
    # Pitcairn Islands
    assert get_gadm_filename("PN") == "gadm41_PCN"
    # Jersey
    assert get_gadm_filename("JE") == "gadm41_JEY"
    # Spratly Islands
    assert get_gadm_filename("XS") == "gadm41_XSP"
    # Guernsey
    assert get_gadm_filename("GG") == "gadm41_GGY"
    # United States Minor Outlying Islands
    assert get_gadm_filename("UM") == "gadm41_UMI"
    # Svalbard islands
    assert get_gadm_filename("SJ") == "gadm41_SJM"
    # Christmas island
    assert get_gadm_filename("CX") == "gadm41_CXR"
    # Afghanistan
    assert get_gadm_filename("AF") == "gadm41_AFG"
    # American Samoa
    assert get_gadm_filename("AS") == "gadm41_ASM"
    # Aruba
    assert get_gadm_filename("AW") == "gadm41_ABW"
    # Germany
    assert get_gadm_filename("DE") == "gadm41_DEU"
    # Micronesia (Federated States of)
    assert get_gadm_filename("FM") == "gadm41_FSM"
    # Micronesia (Federated States of) with different file_prefix
    assert get_gadm_filename("FM", file_prefix="gadm456_") == "gadm456_FSM"


def test_get_gadm_url():
    """
    Verify what is returned by get_gadm_url.
    """
    gadm_filename = get_gadm_filename("FM")
    url_gadm41 = get_gadm_url(
        "https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/",
        gadm_filename,
    )
    assert (
        url_gadm41
        == f"https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/{gadm_filename}.gpkg"
    )


def test_download_gadm():
    """
    Verify what is returned by download_gadm.
    """
    file_prefix_41 = "gadm41_"
    gadm_url_prefix_41 = "https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/"
    gadm_input_file_args_41 = ["data", "gadm"]
    gadm_input_file_gpkg_41, gadm_filename_41 = download_gadm(
        "XK",
        file_prefix_41,
        gadm_url_prefix_41,
        gadm_input_file_args_41,
        update=True,
    )
    assert gadm_input_file_gpkg_41 == get_path(
        path_cwd, "data/gadm/gadm41_XKO/gadm41_XKO.gpkg"
    )
    assert gadm_filename_41 == "gadm41_XKO"
    list_layers_41 = fiona.listlayers(gadm_input_file_gpkg_41)
    assert list_layers_41[0] == "ADM_ADM_0"
    assert list_layers_41[1] == "ADM_ADM_1"
    assert list_layers_41[2] == "ADM_ADM_2"


def test_normed():
    df_input = pd.DataFrame(
        {"A": [1.0, 2.0, 3.0, 4.0, 5.0], "B": [6.0, 7.0, 8.0, 9.0, 10.0]}
    )
    df_output = normed(df_input)
    df_reference = pd.DataFrame(
        {
            "A": [x / 15.0 for x in [1.0, 2.0, 3.0, 4.0, 5.0]],
            "B": [x / 40.0 for x in [6.0, 7.0, 8.0, 9.0, 10.0]],
        }
    )
    df_comparison = df_output.compare(df_reference)
    assert df_comparison.empty
