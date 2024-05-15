# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import pathlib
from scripts._helpers import sets_path_to_root

path_cwd = pathlib.Path.cwd()

def test_sets_path_to_root():
    sets_path_to_root("pypsa-earth")
    current_path = pathlib.Path.cwd()
    assert current_path == path_cwd
