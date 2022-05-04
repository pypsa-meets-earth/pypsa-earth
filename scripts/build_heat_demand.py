# -*- coding: utf-8 -*-
"""Build heat demand time series."""

import os

import atlite
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

if __name__ == "__main__":

    if "snakemake" not in globals():
        from helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_heat_demand", simpl="", clusters="4")
        sets_path_to_root("pypsa-earth-sec")

    time = pd.date_range(freq="h", **snakemake.config["snapshots"])
    cutout_config = snakemake.input.cutout
    cutout = atlite.Cutout(cutout_config).sel(time=time)

    clustered_regions = (
        gpd.read_file(snakemake.input.regions_onshore)
        .set_index("name")
        .buffer(0)
        .squeeze()
    )

    I = cutout.indicatormatrix(clustered_regions)

    for area in ["rural", "urban", "total"]:

        pop_layout = xr.open_dataarray(snakemake.input[f"pop_layout_{area}"])

        stacked_pop = pop_layout.stack(spatial=("y", "x"))
        M = I.T.dot(np.diag(I.dot(stacked_pop)))

        heat_demand = cutout.heat_demand(matrix=M.T, index=clustered_regions.index)

        heat_demand.to_netcdf(snakemake.output[f"heat_demand_{area}"])
