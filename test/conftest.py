# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import pathlib
import shutil

import pandas as pd
import pypsa
import pytest
import yaml

_content_temp_file = "content"
_name_temp_file = "hello.txt"
_temp_content_dir = "temp_content_dir"
_sub_temp_content_dir = "sub_temp_content_dir"


@pytest.fixture(scope="function")
def get_temp_file(tmpdir):
    p = pathlib.Path(tmpdir, _name_temp_file)
    p.write_text(_content_temp_file)
    yield p
    pathlib.Path(p).unlink(missing_ok=True)


@pytest.fixture(scope="function")
def get_temp_folder(tmpdir):
    temp_content_dir = tmpdir.join(_temp_content_dir)
    sub_temp_content_dir = temp_content_dir.join(_sub_temp_content_dir)
    yield sub_temp_content_dir
    shutil.rmtree(str(sub_temp_content_dir))


@pytest.fixture(scope="function")
def get_power_network_scigrid_de():
    return pypsa.examples.scigrid_de()


@pytest.fixture(scope="function")
def get_power_network_ac_dc_meshed():
    return pypsa.examples.ac_dc_meshed()


@pytest.fixture(scope="function")
def get_config_dict():
    path_config = pathlib.Path(pathlib.Path.cwd(), "config.default.yaml")
    with open(path_config, "r") as file:
        config_dict = yaml.safe_load(file)
    return config_dict


@pytest.fixture(scope="function")
def battery_network():
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2022-01-01", periods=24, freq="H"))
    n.add("Bus", "bus1")
    n.add(
        "Link",
        "battery discharger",
        bus0="bus1",
        bus1="bus1",
        p_nom_extendable=True,
        efficiency=0.9,
    )
    n.add("Link", "battery charger", bus0="bus1", bus1="bus1", p_nom_extendable=True)

    # Add links
    n.add(
        "Link",
        "urban central CHP electric",
        bus0="bus0",
        bus1="bus1",
        p_nom_extendable=True,
        efficiency=0.9,
        p_nom_ratio=1.0,
    )
    n.add(
        "Link",
        "urban central CHP heat",
        bus0="bus1",
        bus1="bus0",
        p_nom_extendable=True,
        efficiency=0.9,
    )

    n.optimize.create_model()

    return n
