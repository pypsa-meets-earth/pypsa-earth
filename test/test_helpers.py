# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import pathlib

from scripts._helpers import get_dirname_abs_path, get_path


def test_get_dirname_abs_path():
    """
    Verify the path returned by get_dirname_abs_path()
    """
    path_cwd = str(pathlib.Path.cwd().parent)
    dir_name_file = get_dirname_abs_path(__file__)
    dir_name_cwd = get_dirname_abs_path(".")
    assert dir_name_file == path_cwd + "/test"
    assert dir_name_cwd == path_cwd


def test_get_path():
    """
    Verify the path returned by get_path()
    """
    path_cwd = str(pathlib.Path.cwd().parent)
    path_new = str(get_path(pathlib.Path(__file__).parent, "..", "logs", "rule.log"))
    path = str(get_path(path_cwd, "sub_path_1", "sub_path_2", "sub_path_3", "sub_path_4", "sub_path_5", "file.nc"))
    assert path == path_cwd + "/sub_path_1/sub_path_2/sub_path_3/sub_path_4/sub_path_5/file.nc"
    assert path_new == str(pathlib.Path(__file__).parent.joinpath("..", "logs", "rule.log"))
