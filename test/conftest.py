# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import pathlib

import pytest

_content_temp_file = "content"
_name_temp_file = "hello.txt"


@pytest.fixture(scope="function")
def get_temp_file(tmpdir):
    p = tmpdir.join(_name_temp_file)
    p.write(_content_temp_file)
    yield p
    pathlib.Path(p).unlink(missing_ok=True)
