# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import pathlib
from test.conftest import _content_temp_file, get_temp_file
from scripts._helpers import get_abs_path, get_basename_abs_path, get_dirname_abs_path, get_dirname_path, get_path, get_path_size

path_cwd = str(pathlib.Path.cwd())


def test_get_abs_path():
    """
    Verify the path returned by get_abs_path()
    """
    abs_file = get_abs_path(__file__)
    assert abs_file == __file__

def test_get_dirname_abs_path():
    """
    Verify the path returned by get_dirname_abs_path()
    """
    dir_abs_name_file = get_dirname_abs_path(__file__)
    dir_abs_name_cwd = get_dirname_abs_path(".")
    assert dir_abs_name_file == path_cwd + "/test"
    assert dir_abs_name_cwd == str(pathlib.Path(path_cwd).parent)


def test_get_dirname_path():
    """
    Verify the path returned by get_dirname_path()
    """
    dir_name_file = get_dirname_path(__file__)
    dir_name_cwd = get_dirname_path(".")
    assert dir_name_file == path_cwd + "/test"
    assert dir_name_cwd == "."


def test_get_basename_abs_path():
    """
    Verify the path returned by get_basename_abs_path()
    """
    base_name_file = get_basename_abs_path(__file__)
    assert base_name_file == "test_helpers.py"


def test_get_path():
    """
    Verify the path returned by get_path()
    """
    file_name_path_one = get_path(path_cwd, "sub_path_1", "sub_path_2", "sub_path_3", "sub_path_4", "sub_path_5", "file.nc")
    path_name_path_two = get_path(pathlib.Path(__file__).parent, "..", "logs", "rule.log")
    assert file_name_path_one == path_cwd + "/sub_path_1/sub_path_2/sub_path_3/sub_path_4/sub_path_5/file.nc"
    assert path_name_path_two == str(pathlib.Path(__file__).parent.joinpath("..", "logs", "rule.log"))


def test_get_path_size(get_temp_file):
    """
    Verify the path size (in bytes) returned by get_path_size()
    """
    path = get_temp_file
    file_size = get_path_size(str(path))
    assert file_size == len(_content_temp_file)
