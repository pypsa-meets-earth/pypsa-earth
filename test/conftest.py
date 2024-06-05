# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import pathlib
import shutil

import pytest

_content_temp_file = "content"
_name_temp_file = "hello.txt"
_temp_content_dir = "temp_content_dir"
_sub_temp_content_dir = "sub_temp_content_dir"


@pytest.fixture(scope="function")
def get_temp_file(tmpdir):
    p = tmpdir.join(_name_temp_file)
    p.write(_content_temp_file)
    yield p
    pathlib.Path(p).unlink(missing_ok=True)


@pytest.fixture(scope="function")
def get_temp_folder(tmpdir):
    temp_content_dir = tmpdir.join(_temp_content_dir)
    sub_temp_content_dir = temp_content_dir.join(_sub_temp_content_dir)
    yield sub_temp_content_dir
    shutil.rmtree(str(sub_temp_content_dir))
