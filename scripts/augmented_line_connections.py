# SPDX-FileCopyrightText: : 2021 The PyPSA-Africa Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later
# coding: utf-8
"""
Makes a network more meshed

Relevant Settings
-----------------

.. code:: yaml


.. seealso::


Inputs
------


Outputs
-------



Description
-----------


"""
import logging
import os

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyomo.environ as po
import pypsa
import shapely
from _helpers import configure_logging

logger = logging.getLogger(__name__)


### Functions


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("augmented_line_connections",
                                   network="elec",
                                   simpl="",
                                   clusters="60")
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)
