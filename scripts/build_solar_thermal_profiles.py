# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Build solar thermal collector time series.
"""

import os

import atlite
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_solar_thermal_profiles",
            simpl="",
            clusters="4",
        )

    config = snakemake.params.solar_thermal_config

    time = pd.date_range(freq="h", **snakemake.params.snapshots)
    cutout_config = snakemake.input.cutout
    cutout = atlite.Cutout(cutout_config).sel(time=time)

    clustered_regions = (
        gpd.read_file(snakemake.input.regions_onshore)
        .set_index("name")
        .buffer(0)
        .squeeze()
    )

    I = cutout.indicatormatrix(clustered_regions)

    for area in ["total", "rural", "urban"]:
        pop_layout = xr.open_dataarray(snakemake.input[f"pop_layout_{area}"])

        stacked_pop = pop_layout.stack(spatial=("y", "x"))
        M = I.T.dot(np.diag(I.dot(stacked_pop)))

        nonzero_sum = M.sum(axis=0, keepdims=True)
        nonzero_sum[nonzero_sum == 0.0] = 1.0
        M_tilde = M / nonzero_sum

        solar_thermal = cutout.solar_thermal(
            **config, matrix=M_tilde.T, index=clustered_regions.index
        )

        solar_thermal.to_netcdf(snakemake.output[f"solar_thermal_{area}"])
