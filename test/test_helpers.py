# -*- coding: utf-8 -*-
import os

from scripts._helpers import get_dirname_abs_path


def test_get_dirname_abs_path():
    """
    Verify the path returned by get_dirname_abs_path()
    """
    path_cwd = os.path.dirname(os.getcwd())
    dir_name_file = get_dirname_abs_path(__file__)
    dir_name_cwd = get_dirname_abs_path(".")
    assert dir_name_file == path_cwd + "/test"
    assert dir_name_cwd == path_cwd
