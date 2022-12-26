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
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
from _helpers import mock_snakemake, sets_path_to_root
from build_test_configs import create_test_config
from ruamel.yaml import YAML
from shapely.validation import make_valid


def _multi_index_scen(rulename, keys):
    return pd.MultiIndex.from_product([[rulename], keys], names=["rule", "key"])


def _mock_snakemake(rule, **kwargs):
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    snakemake = mock_snakemake(rule, **kwargs)
    sets_path_to_root("pypsa-earth")
    return snakemake


def generate_scenario_by_country(path_base, country_list):
    "Utility function to generate multiple scenario configs"

    from scripts.download_osm_data import create_country_list

    clean_country_list = create_country_list(country_list)

    for c in clean_country_list:
        create_test_config(path_base, {"countries": [c]}, f"configs/scenarios/{c}.yaml")


def collect_basic_osm_stats(path, rulename, header):
    if os.path.exists(path) and (os.stat(path).st_size > 0):
        df = gpd.read_file(path)
        n_elem = len(df)
    else:
        pd.DataFrame()

    return pd.DataFrame(
        [n_elem], columns=_multi_index_scen(rulename, [header + "-size"])
    )


def collect_network_osm_stats(path, rulename, header, metric_crs="EPSG:3857"):
    if os.path.exists(path) and (os.stat(path).st_size > 0):
        try:
            df = gpd.read_file(path)
            n_elem = len(df)
            obj_length = (
                df["geometry"].apply(make_valid).to_crs(crs=metric_crs).geometry.length
            )
            len_obj = np.nansum(obj_length * df.circuits)

            len_dc_obj = 0.0
            if "frequency" in df.columns:
                coerced_vals = pd.to_numeric(df.frequency, errors="coerce")
                idx_dc = coerced_vals[coerced_vals.as_type(int) == 0].index
                len_dc_obj = obj_length.loc[idx_dc].sum()

            return pd.DataFrame(
                [[n_elem, len_obj, len_dc_obj]],
                columns=_multi_index_scen(
                    rulename,
                    [header + "-" + k for k in ["size", "length", "length_dc"]],
                ),
            )
        except:
            return pd.DataFrame()
    else:
        pd.DataFrame()


def collect_osm_stats(rulename, **kwargs):
    metric_crs = kwargs.pop("metric_crs", "EPSG:3857")
    only_basic = kwargs.pop("only_basic", False)

    df_list = []

    for k, v in kwargs.items():
        if not only_basic and (k in ["lines", "cables"]):
            df_list.append(
                collect_network_osm_stats(v, rulename, k, metric_crs=metric_crs)
            )
        else:
            df_list.append(collect_basic_osm_stats(v, rulename, k))

    return pd.concat(df_list, axis=1)


def collect_raw_osm_stats(rulename="download_osm_data", metric_crs="EPSG:3857"):
    snakemake = _mock_snakemake("download_osm_data")

    options_raw = dict(snakemake.output)
    options_raw.pop("generators_csv")

    return collect_osm_stats(
        rulename, only_basic=True, metric_crs=metric_crs, **options_raw
    )


def collect_clean_osm_stats(rulename="clean_osm_data", metric_crs="EPSG:3857"):
    snakemake = _mock_snakemake("clean_osm_data")

    options_clean = dict(snakemake.output)
    options_clean.pop("generators_csv")

    return collect_osm_stats(rulename, metric_crs=metric_crs, **options_clean)


def collect_network_stats(network_rule, config):

    wildcards = {
        k: str(config["scenario"][k][0]) for k in ["simpl", "clusters", "ll", "opts"]
    }

    snakemake = _mock_snakemake(network_rule, **wildcards)

    network_path = (
        snakemake.output["network"]
        if "network" in snakemake.output.keys()
        else snakemake.output[0]
    )

    if os.path.exists(network_path):
        try:
            n = pypsa.Network(network_path)

            lines_length = (n.lines.length * n.lines.num_parallel).sum()
            lines_capacity = n.lines.s_nom.sum()

            gen_stats = n.generators.groupby("carrier").p_nom.sum()
            hydro_stats = n.storage_units.groupby("carrier").p_nom.sum()

            line_stats = pd.Series(
                [lines_length, lines_capacity], index=["lines_length", "lines_capacity"]
            )

            network_stats = (
                pd.concat([gen_stats, hydro_stats, line_stats]).to_frame().transpose()
            )
            network_stats.columns = _multi_index_scen(
                network_rule, network_stats.columns
            )

            return network_stats
        except:
            return pd.DataFrame()
    else:
        return pd.DataFrame()


def collect_shape_stats(rulename="build_shapes", area_crs="ESRI:54009"):
    snakemake = _mock_snakemake(rulename)

    if not os.path.exists(snakemake.output.africa_shape):
        return pd.DataFrame()

    df_continent = gpd.read_file(snakemake.output.africa_shape)
    continent_area = (
        df_continent["geometry"]
        .apply(make_valid)
        .to_crs(crs=area_crs)
        .geometry.area.iloc[0]
    )

    if not os.path.exists(snakemake.output.gadm_shapes):
        return pd.DataFrame()

    df_gadm = gpd.read_file(snakemake.output.gadm_shapes)
    pop_tot = df_gadm["pop"].sum()
    gdp_tot = df_gadm["gdp"].sum()
    gadm_size = len(df_gadm)
    gadm_country_matching_stats = df_gadm.country.value_counts(normalize=True)
    gadm_country_matching = float(gadm_country_matching_stats.iloc[0]) * 100

    return pd.DataFrame(
        [[continent_area, gadm_size, gadm_country_matching, pop_tot, gdp_tot]],
        columns=_multi_index_scen(
            rulename, ["area", "gadm_size", "country_matching", "pop", "gdp"]
        ),
    )


def collect_snakemake_stats(name, dict_dfs):
    list_rules = [
        "download_osm_data",
        "clean_osm_data",
        "build_shapes",
        "base_network",
        "add_electricity",
        "simplify_network",
        "cluster_network",
        "solve_network",
    ]

    return pd.DataFrame(
        [[rule in dict_dfs.keys() for rule in list_rules]],
        columns=_multi_index_scen(name, list_rules),
    )


def calculate_stats(config, metric_crs="EPSG:3857", area_crs="ESRI:54009"):
    "Function to calculate statistics"

    df_osm_raw = collect_raw_osm_stats(metric_crs=metric_crs)
    df_osm_clean = collect_clean_osm_stats(metric_crs=metric_crs)
    df_shapes = collect_shape_stats(area_crs=area_crs)

    network_dict = {
        network_rule: collect_network_stats(network_rule, config)
        for network_rule in [
            "base_network",
            "add_electricity",
            "simplify_network",
            "cluster_network",
            "solve_network",
        ]
    }

    dict_dfs = {
        "download_osm_data": df_osm_raw,
        "clean_osm_data": df_osm_clean,
        "build_shapes": df_shapes,
        **network_dict,
    }

    df_snakemake = collect_snakemake_stats("snakemake_status", dict_dfs)

    dict_dfs["snakemake_status"] = df_snakemake

    return dict_dfs


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("run_scenario", scenario="NG")

    sets_path_to_root("pypsa-earth")

    # generate_scenario_by_country("configs/scenarios/base.yaml", snakemake.config["countries"])

    scenario = snakemake.wildcards["scenario"]
    dir_scenario = snakemake.output["dir_scenario"]
    stats_scenario = snakemake.output["stats_scenario"]
    base_config = snakemake.config.get("base_config", "./config.default.scenario.yaml")

    # create scenario config
    create_test_config(base_config, f"configs/scenarios/{scenario}.yaml", "config.yaml")

    # execute workflow
    # val = os.system("snakemake -j all solve_all_networks --forceall --rerun-incomplete")

    # create statistics
    stats = calculate_stats(snakemake.config)
    stats = pd.concat(stats.values(), axis=1).set_index(pd.Index([scenario]))
    stats.to_csv(stats_scenario)

    # copy output files
    for f in ["resources", "networks", "results", "benchmarks"]:
        copy_dir = os.path.abspath(f"{dir_scenario}/{f}")
        if os.path.exists(copy_dir):
            shutil.rmtree(copy_dir)
        abs_f = os.path.abspath(f)
        if os.path.exists(abs_f):
            shutil.copytree(abs_f, copy_dir)
            # shutil.rmtree(abs_f, ignore_errors=True)

    shutil.copy("config.yaml", f"{dir_scenario}/config.yaml")
