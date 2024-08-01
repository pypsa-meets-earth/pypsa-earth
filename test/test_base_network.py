# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.append("./scripts")

from _helpers import get_path
from base_network import (
    _get_linetype_by_voltage,
    _get_linetypes_config,
    _load_buses_from_osm,
    _load_converters_from_osm,
    _load_lines_from_osm,
    _load_transformers_from_osm,
    _set_electrical_parameters_converters,
    _set_electrical_parameters_dc_lines,
    _set_electrical_parameters_lines,
    _set_electrical_parameters_links,
    _set_electrical_parameters_transformers,
    get_country,
)

path_cwd = pathlib.Path.cwd()

# Common references

# ---> buses

df_buses_input = pd.DataFrame(
    {
        "bus_id": 0,
        "station_id": 0,
        "voltage": 161000,
        "dc": False,
        "symbol": "substation",
        "under_construction": False,
        "tag_substation": "transmission",
        "tag_area": 0.0,
        "lon": 2.5914,
        "lat": 9.3321,
        "country": "BJ",
        "geometry": "POINT (2.5914 9.3321)",
        "substation_lv": True,
    },
    index=[0],
)

df_buses_reference = pd.DataFrame(
    {
        "bus_id": "0",
        "v_nom": 161.0,
        "symbol": "substation",
        "under_construction": False,
        "tag_substation": "transmission",
        "tag_area": 0.0,
        "lon": 2.5914,
        "lat": 9.3321,
        "country": "BJ",
        "geometry": "POINT (2.5914 9.3321)",
        "substation_lv": True,
        "carrier": "AC",
        "x": 2.5914,
        "y": 9.3321,
    },
    index=[0],
).set_index("bus_id")

# ---> converters

df_converters_input = pd.DataFrame(
    {
        "index": 0,
        "converter_id": "convert_20_41",
        "bus0": "41",
        "bus1": "42",
        "underground": False,
        "under_construction": False,
        "country": "US",
        "geometry": "LINESTRING(-122.3787 37.6821, -122.3777 37.6831)",
    },
    index=[0],
)

df_converters_reference = pd.DataFrame(
    {
        "converter_id": "convert_20_41",
        "Unnamed: 0": 0,
        "index": 0,
        "bus0": "41",
        "bus1": "42",
        "underground": False,
        "under_construction": False,
        "country": "US",
        "geometry": "LINESTRING(-122.3787 37.6821, -122.3777 37.6831)",
        "carrier": "B2B",
        "dc": True,
    },
    index=[0],
).set_index("converter_id")

# ---> lines

df_lines_input = pd.DataFrame(
    {
        "line_id": ["204361221-1_0", "204361287-1_1"],
        "tag_frequency": [50.0, 0.0],
        "tag_type": ["line", "line"],
        "voltage": [161000, 178658],
        "bus0": ["111", "111"],
        "bus1": ["0", "0"],
        "circuits": (3.0, 3.0),
        "length": [110071.89434240988, 118723.89434240988],
        "underground": [False, False],
        "under_construction": [False, False],
        "dc": [False, False],
        "country": ["BJ", "BJ"],
        "geometry": [
            "LINESTRING (2.6594 10.2042, 2.6594451 10.2042341)",
            "LINESTRING (2.6594 10.2042, 2.6594451 10.2042341)",
        ],
        "bounds": [
            "MULTIPOINT ((2.6594 10.2042), (2.5914 9.3321))",
            "MULTIPOINT ((2.6594 10.2042), (2.5914 9.3321))",
        ],
        "bus_0_coors": ["POINT (2.6594 10.2042)", "POINT (2.6594 10.2042)"],
        "bus_1_coors": ["POINT (2.5914 9.3321)", "POINT (2.5914 9.3321)"],
        "bus0_lon": [2.6594, 2.6594],
        "bus0_lat": [10.2042, 10.2042],
        "bus1_lon": [2.5914, 2.5914],
        "bus1_lat": [9.3321, 9.3321],
    }
)

df_lines_reference = pd.DataFrame(
    {
        "line_id": ["204361221-1_0", "204361287-1_1"],
        "tag_frequency": [50.0, 0.0],
        "tag_type": ["line", "line"],
        "v_nom": [161.0, 178.658],
        "bus0": ["111", "111"],
        "bus1": ["0", "0"],
        "num_parallel": [3.0, 3.0],
        "length": [110.07189434240988, 118.72389434240988],
        "underground": [False, False],
        "under_construction": [False, False],
        "dc": [False, False],
        "country": ["BJ", "BJ"],
        "geometry": [
            "LINESTRING (2.6594 10.2042, 2.6594451 10.2042341)",
            "LINESTRING (2.6594 10.2042, 2.6594451 10.2042341)",
        ],
        "bounds": [
            "MULTIPOINT ((2.6594 10.2042), (2.5914 9.3321))",
            "MULTIPOINT ((2.6594 10.2042), (2.5914 9.3321))",
        ],
        "bus_0_coors": ["POINT (2.6594 10.2042)", "POINT (2.6594 10.2042)"],
        "bus_1_coors": ["POINT (2.5914 9.3321)", "POINT (2.5914 9.3321)"],
        "bus0_lon": [2.6594, 2.6594],
        "bus0_lat": [10.2042, 10.2042],
        "bus1_lon": [2.5914, 2.5914],
        "bus1_lat": [9.3321, 9.3321],
    }
).set_index("line_id")

lines_ac_reference = pd.DataFrame(
    {
        "tag_frequency": 50.0,
        "tag_type": "line",
        "v_nom": 161.0,
        "bus0": "111",
        "bus1": "0",
        "num_parallel": 3.0,
        "length": 110.07189434240988,
        "underground": False,
        "under_construction": False,
        "dc": False,
        "country": "BJ",
        "geometry": "LINESTRING (2.6594 10.2042, 2.6594451 10.2042341)",
        "bounds": "MULTIPOINT ((2.6594 10.2042), (2.5914 9.3321))",
        "bus_0_coors": "POINT (2.6594 10.2042)",
        "bus_1_coors": "POINT (2.5914 9.3321)",
        "bus0_lon": 2.6594,
        "bus0_lat": 10.2042,
        "bus1_lon": 2.5914,
        "bus1_lat": 9.3321,
        "carrier": "AC",
        "type": "243-AL1/39-ST1A 20.0",
        "s_max_pu": 0.7,
    },
    index=[0],
).set_index("tag_frequency")

lines_dc_reference = pd.DataFrame(
    {
        "tag_frequency": 0.0,
        "tag_type": "line",
        "v_nom": 178.658,
        "bus0": "111",
        "bus1": "0",
        "num_parallel": 3.0,
        "length": 118.72389434240988,
        "underground": False,
        "under_construction": False,
        "dc": True,
        "country": "BJ",
        "geometry": "LINESTRING (2.6594 10.2042, 2.6594451 10.2042341)",
        "bounds": "MULTIPOINT ((2.6594 10.2042), (2.5914 9.3321))",
        "bus_0_coors": "POINT (2.6594 10.2042)",
        "bus_1_coors": "POINT (2.5914 9.3321)",
        "bus0_lon": 2.6594,
        "bus0_lat": 10.2042,
        "bus1_lon": 2.5914,
        "bus1_lat": 9.3321,
        "carrier": "DC",
        "type": "HVDC XLPE 1000",
        "s_max_pu": 0.7,
    },
    index=[0],
).set_index("tag_frequency")

lines_dict = {
    "ac_types": {
        132.0: "243-AL1/39-ST1A 20.0",
        220.0: "Al/St 240/40 2-bundle 220.0",
        300.0: "Al/St 240/40 3-bundle 300.0",
        380.0: "Al/St 240/40 4-bundle 380.0",
        500.0: "Al/St 240/40 4-bundle 380.0",
        750.0: "Al/St 560/50 4-bundle 750.0",
    },
    "dc_types": {
        500.0: "HVDC XLPE 1000",
    },
    "s_max_pu": 0.7,
    "s_nom_max": np.inf,
    "length_factor": 1.25,
    "under_construction": "zero",
}

# ---> links

links_dict = {
    "p_max_pu": 2.1,
    "p_nom_max": np.inf,
    "under_construction": "zero",
}

# ---> transformers

transformers_dict = {
    "x": 0.1,
    "s_nom": 2000.0,
    "type": "",
}

df_transformers_input = pd.DataFrame(
    {
        "line_id": "transf_1_0",
        "bus0": "1",
        "bus1": "2",
        "voltage_bus0": 161000,
        "voltage_bus1": 330000,
        "country": "BJ",
        "geometry": "LINESTRING(2.648 6.7394, 2.649 6.7404)",
        "bounds": "MULTIPOINT((2.648 6.7394), (2.649 6.7404))",
        "bus_0_coors": "POINT(2.648 6.7394)",
        "bus_1_coors": "POINT(2.649 6.7404)",
        "bus0_lon": 2.648,
        "bus0_lat": 6.7394,
        "bus1_lon": 2.649,
        "bus1_lat": 6.7404,
    },
    index=[0],
)

df_transformers_reference = pd.DataFrame(
    {
        "transformer_id": "transf_1_0",
        "Unnamed: 0": 0,
        "bus0": "1",
        "bus1": "2",
        "voltage_bus0": 161000,
        "voltage_bus1": 330000,
        "country": "BJ",
        "geometry": "LINESTRING(2.648 6.7394, 2.649 6.7404)",
        "bounds": "MULTIPOINT((2.648 6.7394), (2.649 6.7404))",
        "bus_0_coors": "POINT(2.648 6.7394)",
        "bus_1_coors": "POINT(2.649 6.7404)",
        "bus0_lon": 2.648,
        "bus0_lat": 6.7394,
        "bus1_lon": 2.649,
        "bus1_lat": 6.7404,
    },
    index=[0],
).set_index("transformer_id")

# ---> voltages

voltages_list = [132.0, 220.0, 300.0, 380.0, 500.0, 750.0]


def test_get_country():
    """
    Verify what returned by get_country()
    """
    data_list = [['"country"=>"NG"'], ['"country"=>"CH"'], ['"country"=>"AU"']]
    df_exercise_with_tags = pd.DataFrame(data_list, columns=["tags"])
    df_exercise_no_tags = pd.DataFrame(data_list, columns=["other"])
    series_with_tags = get_country(df_exercise_with_tags)
    reference_series_with_tags = pd.Series(["NG", "CH", "AU"])
    comparison_series_with_tags = series_with_tags.compare(reference_series_with_tags)
    series_no_tags = get_country(df_exercise_no_tags)
    reference_series_no_tags = pd.Series([np.nan, np.nan, np.nan])
    comparison_series_no_tags = series_no_tags.compare(reference_series_no_tags)
    assert comparison_series_with_tags.size == 0
    assert comparison_series_no_tags.size == 0


def test_load_buses_from_osm(tmpdir):
    """
    Verify what returned by _load_buses_from_osm.
    """
    file_path = get_path(tmpdir, "buses_exercise.csv")
    df_buses_input.to_csv(file_path)
    df_buses_output = _load_buses_from_osm(file_path)
    df_buses_comparison = df_buses_output.compare(df_buses_reference)
    pathlib.Path.unlink(file_path)
    assert df_buses_comparison.empty


def test_load_lines_from_osm(tmpdir):
    """
    Verify what returned by _load_lines_from_osm.
    """
    file_path = get_path(tmpdir, "lines_exercise.csv")
    df_lines_input.to_csv(file_path)
    df_lines_output = _load_lines_from_osm(file_path)
    df_lines_comparison = df_lines_output.compare(df_lines_reference)
    pathlib.Path.unlink(file_path)
    assert df_lines_comparison.empty


def test_load_transformers_from_osm(tmpdir):
    """
    Verify what returned by _load_transformers_from_osm.
    """
    file_path = get_path(tmpdir, "transformers_exercise.csv")
    df_transformers_input.to_csv(file_path)
    df_transformers_output = _load_transformers_from_osm(file_path)
    df_transformers_comparison = df_transformers_output.compare(
        df_transformers_reference
    )
    pathlib.Path.unlink(file_path)
    assert df_transformers_comparison.empty


def test_load_converters_from_osm(tmpdir):
    """
    Verify what returned by _load_converters_from_osm.
    """
    file_path = get_path(tmpdir, "converters_exercise.csv")
    df_converters_input.to_csv(file_path)
    df_converters_output = _load_converters_from_osm(file_path)
    df_converters_comparison = df_converters_output.compare(df_converters_reference)
    pathlib.Path.unlink(file_path)
    assert df_converters_comparison.empty


def test_get_linetypes_config():
    """
    Verify what returned by _get_linetypes_config.
    """
    output_dict_ac = _get_linetypes_config(lines_dict["ac_types"], voltages_list)
    output_dict_dc = _get_linetypes_config(lines_dict["dc_types"], voltages_list)
    assert output_dict_ac == lines_dict["ac_types"]
    assert output_dict_dc == lines_dict["dc_types"]


def test_get_linetype_by_voltage():
    """
    Verify what returned by _get_linetype_by_voltage.
    """
    v_nom_list = [
        50.0,
        101.0,
        180.0,
        210.0,
        220.0,
        225.0,
        285.0,
        300.0,
        333.0,
        390.0,
        600.0,
        750.0,
        800.0,
    ]

    line_type_list = []

    for v_nom in v_nom_list:
        line_type_list.append(_get_linetype_by_voltage(v_nom, lines_dict["ac_types"]))

    assert line_type_list == [
        "243-AL1/39-ST1A 20.0",
        "243-AL1/39-ST1A 20.0",
        "Al/St 240/40 2-bundle 220.0",
        "Al/St 240/40 2-bundle 220.0",
        "Al/St 240/40 2-bundle 220.0",
        "Al/St 240/40 2-bundle 220.0",
        "Al/St 240/40 3-bundle 300.0",
        "Al/St 240/40 3-bundle 300.0",
        "Al/St 240/40 3-bundle 300.0",
        "Al/St 240/40 4-bundle 380.0",
        "Al/St 240/40 4-bundle 380.0",
        "Al/St 560/50 4-bundle 750.0",
        "Al/St 560/50 4-bundle 750.0",
    ]


def test_set_electrical_parameters_lines(tmpdir):
    """
    Verify what returned by _set_electrical_parameters_lines.
    """
    file_path = get_path(tmpdir, "lines_exercise.csv")
    df_lines_input.to_csv(file_path)
    df_lines_output = _load_lines_from_osm(file_path).reset_index(drop=True)
    df_lines_output_ac = df_lines_output[
        df_lines_output.tag_frequency.astype(float) != 0
    ].copy()
    df_lines_output_dc = df_lines_output[
        df_lines_output.tag_frequency.astype(float) == 0
    ].copy()
    lines_ac = _set_electrical_parameters_lines(
        lines_dict, voltages_list, df_lines_output_ac
    ).set_index("tag_frequency")
    lines_dc = _set_electrical_parameters_dc_lines(
        lines_dict, voltages_list, df_lines_output_dc
    ).set_index("tag_frequency")
    df_lines_ac_comparison = lines_ac.compare(lines_ac_reference)
    df_lines_dc_comparison = lines_dc.compare(lines_dc_reference)
    pathlib.Path.unlink(file_path)
    assert df_lines_ac_comparison.empty
    assert df_lines_dc_comparison.empty


def test_set_electrical_parameters_links(tmpdir):
    """
    Verify what returned by _set_electrical_parameters_links.
    """
    file_path = get_path(tmpdir, "lines_exercise.csv")
    df_lines_input.to_csv(file_path)
    df_lines_output = _load_lines_from_osm(file_path).reset_index(drop=True)
    df_lines_output_dc = df_lines_output[
        df_lines_output.tag_frequency.astype(float) == 0
    ].copy()
    lines_dc = _set_electrical_parameters_dc_lines(
        lines_dict, voltages_list, df_lines_output_dc
    )
    new_lines_dc = _set_electrical_parameters_links(links_dict, lines_dc).set_index(
        "tag_frequency"
    )
    new_lines_dc_reference = lines_dc_reference.copy(deep=True)
    new_lines_dc_reference["p_max_pu"] = links_dict["p_max_pu"]
    new_lines_dc_reference["p_min_pu"] = -links_dict["p_max_pu"]
    pathlib.Path.unlink(file_path)
    df_comparison = new_lines_dc.compare(new_lines_dc_reference)
    assert df_comparison.empty


def test_set_electrical_parameters_transformers(tmpdir):
    """
    Verify what returned by _set_electrical_parameters_transformers.
    """
    file_path = get_path(tmpdir, "transformers_exercise.csv")
    df_transformers_input.to_csv(file_path)
    df_transformers_output = _load_transformers_from_osm(file_path)
    df_transformers_parameters = _set_electrical_parameters_transformers(
        transformers_dict, df_transformers_output
    )
    df_transformers_parameters_reference = df_transformers_reference.copy(deep=True)
    df_transformers_parameters_reference["x"] = transformers_dict["x"]
    df_transformers_parameters_reference["s_nom"] = transformers_dict["s_nom"]
    df_transformers_parameters_reference["type"] = transformers_dict["type"]
    pathlib.Path.unlink(file_path)
    df_comparison = df_transformers_parameters.compare(
        df_transformers_parameters_reference
    )
    assert df_comparison.empty


def test_set_electrical_parameters_converters(tmpdir):
    """
    Verify what returned by _set_electrical_parameters_converters.
    """
    file_path = get_path(tmpdir, "converters_exercise.csv")
    df_converters_input.to_csv(file_path)
    df_converters_output = _load_converters_from_osm(file_path)
    df_converters_parameters = _set_electrical_parameters_converters(
        links_dict, df_converters_output
    )
    df_converters_parameters_reference = df_converters_reference.copy(deep=True)
    df_converters_parameters_reference["p_max_pu"] = links_dict["p_max_pu"]
    df_converters_parameters_reference["p_min_pu"] = -links_dict["p_max_pu"]
    df_converters_parameters_reference["p_nom"] = 2000
    df_converters_parameters_reference["under_construction"] = False
    df_converters_parameters_reference["underground"] = False
    pathlib.Path.unlink(file_path)
    df_comparison = df_converters_parameters.compare(df_converters_parameters_reference)
    assert df_comparison.empty
