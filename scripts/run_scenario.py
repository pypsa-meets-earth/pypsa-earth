# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2022 PyPSA-Earth Authors
#
# coding: utf-8
"""
Execute a scenario optimization

Run iteratively the workflow under different conditions
and store the results in specific folders
"""
import collections.abc
import copy
import os
import shutil

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
from _helpers import sets_path_to_root
from build_test_configs import create_test_config
from shapely.validation import make_valid


def generate_scenario_by_country(path_base, country_list):
    "Utility function to generate multiple scenario configs"

    from scripts.download_osm_data import create_country_list

    clean_country_list = create_country_list(country_list)

    for c in clean_country_list:
        create_test_config(path_base, {"countries": [c]}, f"configs/scenarios/{c}.yaml")


def collect_basic_osm_stats(path, name):
    if os.path.exists(path):
        df = gpd.read_file(path)
        n_elem = len(df)
    else:
        n_elem = np.nan

    return pd.DataFrame({"size": n_elem}, index=[name])


def collect_network_osm_stats(path, name, crs_metric=3857):
    if os.path.exists(path):
        df = gpd.read_file(path)
        n_elem = len(df)
        obj_length = (
            df["geometry"].apply(make_valid).to_crs(epsg=crs_metric).geometry.length
        )
        len_obj = sum(obj_length)

        len_dc_obj = 0.0
        if "frequency" in df.columns:
            coerced_vals = pd.to_numeric(df.frequency, errors="coerce")
            idx_dc = coerced_vals[coerced_vals.as_type(int) == 0].index
            len_dc_obj = obj_length.loc[idx_dc].sum()
    else:
        n_elem = np.nan
        len_obj = np.nan
        len_dc_obj = np.nan

    return pd.DataFrame(
        {"size": n_elem, "length": len_obj, "length_dc": len_dc_obj}, index=[name]
    )


def collect_osm_stats(**kwargs):
    crs_metric = kwargs.pop("crs_metric", 3857)
    only_basic = kwargs.pop("only_basic", False)

    df_list = []

    for k, v in kwargs.items():
        if not only_basic and (k in ["lines", "cables"]):
            df_list.append(collect_network_osm_stats(v, k, crs_metric=crs_metric))
        else:
            df_list.append(collect_basic_osm_stats(v, k))

    return pd.concat(df_list)


def collect_raw_osm_stats(crs_metric=3857):
    from _helpers import mock_snakemake

    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    snakemake = mock_snakemake("download_osm_data")

    options_raw = dict(snakemake.output)
    options_raw.pop("generators_csv")

    return collect_osm_stats(only_basic=True, crs_metric=crs_metric, **options_raw)


def collect_clean_osm_stats(crs_metric=3857):
    from _helpers import mock_snakemake

    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    snakemake = mock_snakemake("clean_osm_data")

    options_clean = dict(snakemake.output)
    options_clean.pop("generators_csv")

    return collect_osm_stats(crs_metric=crs_metric, **options_clean)


def collect_network_stats(n_path):
    if os.path.exists(n_path):
        n = pypsa.Network(n_path)
        return n.statistics()
    else:
        return pd.DataFrame()


def collect_shape_stats(area_crs="ESRI:54009"):
    from _helpers import mock_snakemake

    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    snakemake = mock_snakemake("build_shapes")

    continent_area = np.nan
    if os.path.exists(snakemake.output.africa_shape):
        df_continent = gpd.read_file(snakemake.output.africa_shape)
        continent_area = (
            df_continent["geometry"]
            .apply(make_valid)
            .to_crs(crs=area_crs)
            .geometry.area.iloc[0]
        )

    pop_tot = np.nan
    gdp_tot = np.nan
    gadm_size = np.nan
    if os.path.exists(snakemake.output.gadm_shapes):
        df_gadm = gpd.read_file(snakemake.output.gadm_shapes)
        pop_tot = df_gadm["pop"].sum()
        gdp_tot = df_gadm["gdp"].sum()
        gadm_size = len(df_gadm)

    return pd.DataFrame(
        {
            "area": [continent_area],
            "gadm_size": [gadm_size],
            "pop": [pop_tot],
            "gdp": [gdp_tot],
        }
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("run_scenario", scenario="NG")

    sets_path_to_root("pypsa-earth")

    collect_shape_stats()

    # generate_scenario_by_country("configs/scenarios/base.yaml", snakemake.config["countries"])

    scenario = snakemake.wildcards["scenario"]
    base_config = snakemake.config.get("base_config", "./config.default.yaml")

    # create scenario config
    create_test_config(base_config, f"configs/scenarios/{scenario}.yaml", "config.yaml")

    val = os.system("snakemake -j all solve_all_networks --forceall")

    for f in ["resources", "networks", "results", "benchmarks"]:
        shutil.copytree(f, f"scenarios/{scenario}/{f}")

    shutil.copy("config.yaml", f"scenarios/{scenario}/config.yaml")
