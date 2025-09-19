# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: PyPSA-ASEAN, PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Append the CAGR rate estimation using data sourced from https://www.aseanstats.org/
"""
import pandas as pd
from _helpers import read_csv_nafix

efficiency_gains_cagr = read_csv_nafix(
    snakemake.input[0], index_col=0
)
