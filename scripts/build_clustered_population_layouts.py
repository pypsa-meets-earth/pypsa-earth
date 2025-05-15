# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Build clustered population layouts.
"""
import os

import atlite
import geopandas as gpd
import pandas as pd
import xarray as xr
from _helpers import read_csv_nafix, to_csv_nafix

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_clustered_population_layouts",
            simpl="",
            clusters=4,
        )

    cutout_path = (
        snakemake.input.cutout
    )  # os.path.abspath(snakemake.config["atlite"]["cutout"])
    cutout = atlite.Cutout(cutout_path)
    # cutout = atlite.Cutout(snakemake.config['atlite']['cutout'])

    clustered_regions = (
        gpd.read_file(snakemake.input.regions_onshore)
        .set_index("name")
        .buffer(0)
        .squeeze()
    )

    I = cutout.indicatormatrix(clustered_regions)

    pop = {}
    for item in ["total", "urban", "rural"]:
        pop_layout = xr.open_dataarray(snakemake.input[f"pop_layout_{item}"])
        pop[item] = I.dot(pop_layout.stack(spatial=("y", "x")))

    pop = pd.DataFrame(pop, index=clustered_regions.index)

    pop["ct"] = gpd.read_file(snakemake.input.regions_onshore).set_index("name").country
    country_population = pop.total.groupby(pop.ct).sum()
    pop["fraction"] = (pop.total / pop.ct.map(country_population)).fillna(0.0)

    to_csv_nafix(pop, snakemake.output.clustered_pop_layout)

    gdp_layout = xr.open_dataarray(snakemake.input["gdp_layout"])
    gdp = I.dot(gdp_layout.stack(spatial=("y", "x")))
    gdp = pd.DataFrame(gdp, index=clustered_regions.index, columns=["total"])

    gdp["ct"] = gpd.read_file(snakemake.input.regions_onshore).set_index("name").country
    country_gdp = gdp.total.groupby(gdp.ct).sum()
    gdp["fraction"] = (gdp.total / gdp.ct.map(country_gdp)).fillna(0.0)
    to_csv_nafix(gdp, snakemake.output.clustered_gdp_layout)
