# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Create statistics for a given scenario run.

This script contains functions to create statistics of the workflow for the current execution

Relevant statistics that are created are:

- For clean_osm_data and download_osm_data,
  the number of elements, length of the lines and length of dc lines are stored
- For build_shapes, the surface, total GDP, total population and number of shapes are collected
- For build_renewable_profiles, total available potential and average production are collected
- For network rules (base_network, add_electricity, simplify_network and solve_network),
  length of lines, number of buses and total installed capacity by generation technology
- Execution time for the rules, when benchmark is available

Outputs
-------
This rule creates a dataframe containing in the columns the relevant statistics for the current run.
"""
import os
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import create_logger, mock_snakemake, to_csv_nafix
from build_test_configs import create_test_config
from shapely.validation import make_valid

logger = create_logger(__name__)


def _multi_index_scen(rulename, keys):
    return pd.MultiIndex.from_product([[rulename], keys], names=["rule", "key"])


def _mock_snakemake(rule, **kwargs):

    snakemake = mock_snakemake(rule, **kwargs)

    return snakemake


def generate_scenario_by_country(
    path_base, country_list, out_dir="configs/scenarios", pre="config."
):
    """
    Utility function to create copies of a standard yaml file available in
    path_base for every country in country_list. Copies are saved into the
    output directory out_dir.

    Note:
    - the clusters are automatically modified for selected countries with limited data
    - for landlocked countries, offwind technologies are removed (solar, onwind and hydro are forced)

    Parameters
    ----------
        path_base : str
            Path to the standard yaml file used as default
        country_list : list
            List of countries.
            Note: the input is parsed using download_osm_data.create_country_list
        out_dir : str (optional)
            Output directory where output configuration files are executed
    """
    from _helpers import create_country_list, three_2_two_digits_country

    clean_country_list = create_country_list(country_list)

    # file available from https://worldpopulationreview.com/country-rankings/landlocked-countries
    df_landlocked = pd.read_csv("landlocked.csv")
    df_landlocked["countries"] = df_landlocked.cca2.map(three_2_two_digits_country)

    n_clusters = {
        "MG": 3,  # Africa
        "BF": 1,
        "BI": 1,
        "BJ": 2,
        "DJ": 1,
        "GM": 2,
        "LR": 2,
        "LS": 3,
        "NE": 4,
        "SL": 1,
        "SZ": 4,
        "TG": 1,
        "CG": 2,
        "GN": 3,
        "SS": 1,
        "ML": 1,
        "TD": 2,
        "CF": 1,
        "ER": 1,  # Africa
        "GF": 3,  # South America
        "SR": 1,  # South America
        "SG": 1,  # Asia
        "FJ": 4,  # Oceania
    }

    for c in clean_country_list:
        modify_dict = {"countries": [c]}

        if c in n_clusters.keys():
            modify_dict["scenario"] = {"clusters": [str(n_clusters[c])]}

        if df_landlocked["countries"].str.contains(c).any():
            modify_dict["electricity"] = {
                "renewable_carriers": ["solar", "onwind", "hydro"]
            }

        create_test_config(path_base, modify_dict, f"{out_dir}/{pre}{c}.yaml")


def collect_basic_osm_stats(path, rulename, header):
    """
    Collect basic statistics on OSM data: number of items
    """
    if Path(path).is_file() and Path(path).stat().st_size > 0:
        df = gpd.read_file(path)
        n_elem = len(df)

        return pd.DataFrame(
            [n_elem], columns=_multi_index_scen(rulename, [header + "-size"])
        )

    else:
        return pd.DataFrame()


def collect_network_osm_stats(path, rulename, header, metric_crs="EPSG:3857"):
    """
    Collect statistics on OSM network data:
    - number of items
    - length of the stored shapes
    - length of objects with tag_frequency == 0 (DC elements)
    """
    if Path(path).is_file() and Path(path).stat().st_size > 0:
        df = gpd.read_file(path)
        n_elem = len(df)
        obj_length = (
            df["geometry"].apply(make_valid).to_crs(crs=metric_crs).geometry.length
        )
        len_obj = np.nansum(obj_length * df.circuits)

        len_dc_obj = 0.0
        n_elem_dc = 0
        if "tag_frequency" in df.columns:
            coerced_vals = pd.to_numeric(df.tag_frequency, errors="coerce")
            coerced_vals.fillna(50, inplace=True)
            idx_dc = coerced_vals[coerced_vals.astype(int) == 0].index
            len_dc_obj = np.nansum(obj_length.loc[idx_dc] * df.circuits.loc[idx_dc])
            n_elem_dc = int(len(idx_dc))

        return pd.DataFrame(
            [[n_elem, n_elem_dc, len_obj, len_dc_obj]],
            columns=_multi_index_scen(
                rulename,
                [header + "-" + k for k in ["size", "size_dc", "length", "length_dc"]],
            ),
        )
    else:
        return pd.DataFrame()


def collect_osm_stats(rulename, **kwargs):
    """
    Collect statistics on OSM data.

    When lines and cables are considered, then network-related
    statistics are collected (collect_network_osm_stats), otherwise
    basic statistics are (collect_basic_osm_stats)
    """
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
    """
    Collect basic statistics on OSM data; used for raw OSM data.
    """
    snakemake = _mock_snakemake("download_osm_data")

    options_raw = dict(snakemake.output)
    options_raw.pop("generators_csv")

    df_raw_osm_stats = collect_osm_stats(
        rulename, only_basic=True, metric_crs=metric_crs, **options_raw
    )

    add_computational_stats(df_raw_osm_stats, snakemake)

    return df_raw_osm_stats


def collect_clean_osm_stats(rulename="clean_osm_data", metric_crs="EPSG:3857"):
    """
    Collect statistics on OSM data; used for clean OSM data.
    """
    snakemake = _mock_snakemake("clean_osm_data")

    options_clean = dict(snakemake.output)
    options_clean.pop("generators_csv")

    df_clean_osm_stats = collect_osm_stats(
        rulename, metric_crs=metric_crs, **options_clean
    )

    add_computational_stats(df_clean_osm_stats, snakemake)

    return df_clean_osm_stats


def collect_bus_regions_stats(bus_region_rule="build_bus_regions"):
    """
    Collect statistics on bus regions.

    - number of onshore regions
    - number of offshore regions
    """
    snakemake = _mock_snakemake(bus_region_rule)

    fp_onshore = snakemake.output.regions_onshore
    fp_offshore = snakemake.output.regions_offshore

    df = pd.DataFrame()

    if Path(fp_onshore).is_file() and Path(fp_offshore).is_file():
        gdf_onshore = gpd.read_file(fp_onshore)
        gdf_offshore = gpd.read_file(fp_offshore)

        df = pd.DataFrame(
            [[gdf_onshore.shape[0], gdf_offshore.shape[0]]],
            columns=_multi_index_scen(
                bus_region_rule,
                ["n_onshore", "n_offshore"],
            ),
        )

    add_computational_stats(df, snakemake)

    return df


def collect_network_stats(network_rule, scenario_config):
    """
    Collect statistics on pypsa networks:
    - installed capacity by carrier
    - lines total length (accounting for parallel lines)
    - lines total capacity
    """
    wildcards = {
        k: str(scenario_config[k][0]) for k in ["simpl", "clusters", "ll", "opts"]
    }

    snakemake = _mock_snakemake(network_rule, **wildcards)

    network_path = (
        snakemake.output["network"]
        if "network" in snakemake.output.keys()
        else snakemake.output[0]
    )

    def capacity_stats(df):
        if df.empty:
            return pd.Series(dtype=float)
        else:
            return df.groupby("carrier").p_nom.sum().astype(float)

    if Path(network_path).is_file():
        n = pypsa.Network(network_path)

        lines_length = float((n.lines.length * n.lines.num_parallel).sum())

        lines_capacity = float(n.lines.s_nom.sum())
        buses_number = int(n.buses.shape[0])
        buses_number_dc = (n.buses.carrier != "AC").sum()

        network_stats = pd.DataFrame(
            [[buses_number, buses_number_dc, lines_length, lines_capacity]],
            columns=_multi_index_scen(
                network_rule,
                [
                    "buses_number",
                    "buses_number_non_ac",
                    "lines_length",
                    "lines_capacity",
                ],
            ),
        )

        gen_stats = pd.concat(
            [capacity_stats(n.generators), capacity_stats(n.storage_units)], axis=0
        )

        # add demand to statistics
        tot_demand = float(n.loads_t.p_set.sum(axis=1).mean()) * 8760
        df_demand = pd.DataFrame({"demand": [tot_demand]})
        df_demand.columns = _multi_index_scen(network_rule, df_demand.columns)
        network_stats = pd.concat([network_stats, df_demand], axis=1)

        if not gen_stats.empty:
            df_gen_stats = gen_stats.to_frame().transpose().reset_index(drop=True)
            df_gen_stats.columns = _multi_index_scen(network_rule, df_gen_stats.columns)
            network_stats = pd.concat([network_stats, df_gen_stats], axis=1)

        add_computational_stats(network_stats, snakemake)

        return network_stats
    else:
        return pd.DataFrame()


def collect_shape_stats(rulename="build_shapes", area_crs="ESRI:54009"):
    """
    Collect statistics on the shapes created by the workflow:
    - area
    - number of gadm shapes
    - Percentage of shapes having country flag matching the gadm file
    - total population
    - total gdp
    """
    snakemake = _mock_snakemake(rulename)

    if not Path(snakemake.output.africa_shape).is_file():
        return pd.DataFrame()

    df_continent = gpd.read_file(snakemake.output.africa_shape)
    continent_area = (
        df_continent["geometry"]
        .apply(make_valid)
        .to_crs(crs=area_crs)
        .geometry.area.iloc[0]
    )

    if not Path(snakemake.output.gadm_shapes).is_file():
        return pd.DataFrame()

    df_gadm = gpd.read_file(snakemake.output.gadm_shapes)
    pop_tot = float(df_gadm["pop"].sum())
    gdp_tot = float(df_gadm["gdp"].sum())
    gadm_size = len(df_gadm)
    gadm_country_matching_stats = df_gadm.country.value_counts(normalize=True)
    gadm_country_matching = float(gadm_country_matching_stats.iloc[0]) * 100

    df_shapes_stats = pd.DataFrame(
        [[continent_area, gadm_size, gadm_country_matching, pop_tot, gdp_tot]],
        columns=_multi_index_scen(
            rulename, ["area", "gadm_size", "country_matching", "pop", "gdp"]
        ),
    )

    add_computational_stats(df_shapes_stats, snakemake)

    return df_shapes_stats


def collect_only_computational(rulename):
    """
    Rule to create only computational statistics of rule rulename.
    """
    snakemake = _mock_snakemake(rulename)
    df = pd.DataFrame()
    add_computational_stats(df, snakemake)
    return df


def collect_snakemake_stats(
    name, dict_dfs, renewable_config, renewable_carriers_config
):
    """
    Collect statistics on what rules have been successful.
    """
    ren_techs = [tech for tech in renewable_config if tech in renewable_carriers_config]

    list_rules = [
        "download_osm_data",
        "clean_osm_data",
        "build_shapes",
        "build_osm_network",
        "build_bus_regions",
        "build_demand_profiles",
        "build_powerplants",
        *[f"build_renewable_profiles_{rtech}" for rtech in ren_techs],
        "base_network",
        "add_electricity",
        "simplify_network",
        "cluster_network",
        "solve_network",
    ]

    return pd.DataFrame(
        [
            [
                (rule in dict_dfs.keys()) and not dict_dfs[rule].empty
                for rule in list_rules
            ]
        ],
        columns=_multi_index_scen(name, list_rules),
    )


def aggregate_computational_stats(name, dict_dfs):
    """
    Function to aggregate the total computational statistics of the rules.
    """
    cols_comp = ["total_time", "mean_load", "max_memory"]

    def get_selected_cols(df, level=1, lvl_cols=cols_comp):
        if df.columns.nlevels < level + 1:
            return pd.DataFrame()

        filter_cols = df.columns[df.columns.get_level_values(level).isin(lvl_cols)]

        if len(filter_cols) == len(lvl_cols):
            df_ret = df[filter_cols].copy()
            df_ret.columns = filter_cols.get_level_values(level)
            return df_ret
        else:
            return pd.DataFrame()

    def weigh_avg(df, coldata="total_time", colweight="mean_load"):
        d, w = df[coldata], df[colweight]
        return (d * w).sum() / w.sum()

    df_comb = pd.concat(
        [get_selected_cols(d) for d in dict_dfs.values()], ignore_index=True
    )

    if df_comb.empty:
        return pd.DataFrame()

    df_comb_agg = df_comb.agg({"total_time": np.sum, "max_memory": np.max})
    df_comb_agg["mean_load"] = weigh_avg(df_comb)

    df_comb_agg = df_comb_agg.to_frame().transpose().reset_index(drop=True)
    df_comb_agg.columns = _multi_index_scen(name, df_comb_agg.columns)

    return df_comb_agg


def collect_renewable_stats(rulename, technology):
    """
    Collect statistics on the renewable time series generated by the workflow:
    - potential
    - average production by plant (hydro) or bus (other RES)
    """
    snakemake = _mock_snakemake(rulename, technology=technology)

    if Path(snakemake.output.profile).is_file():
        res = xr.open_dataset(snakemake.output.profile)

        if technology == "hydro":
            potential = float(res.inflow.sum())
            avg_production_pu = float(res.inflow.mean(dim=["plant"]).sum())
        else:
            potential = float(res.potential.sum())
            avg_production_pu = float(res.profile.mean(dim=["bus"]).sum())

        df_RES_stats = pd.DataFrame(
            [[potential, avg_production_pu]],
            columns=_multi_index_scen(
                f"{rulename}_{technology}",
                ["potential", "avg_production_pu"],
            ),
        )

        add_computational_stats(df_RES_stats, snakemake, f"{rulename}_{technology}")

        return df_RES_stats
    else:
        return pd.DataFrame()


def add_computational_stats(df, snakemake, column_name=None):
    """
    Add the major computational information of a given rule into the existing
    dataframe.
    """
    comp_data = [np.nan] * 3  # total_time, mean_load and max_memory

    if snakemake.benchmark:
        if not Path(snakemake.benchmark).is_file():
            return df

        bench_data = pd.read_csv(snakemake.benchmark, delimiter="\t")

        comp_data = bench_data[["s", "mean_load", "max_vms"]].iloc[0].values

    if column_name is None:
        column_name = snakemake.rule

    new_cols = _multi_index_scen(column_name, ["total_time", "mean_load", "max_memory"])

    df[new_cols] = [comp_data]

    if len(df.columns) == len(new_cols):
        df.columns = new_cols

    return df


def calculate_stats(
    scenario_config,
    renewable_config,
    renewable_carriers_config,
    metric_crs="EPSG:3857",
    area_crs="ESRI:54009",
):
    "Function to collect all statistics"
    df_osm_raw = collect_raw_osm_stats(metric_crs=metric_crs)
    df_osm_clean = collect_clean_osm_stats(metric_crs=metric_crs)
    df_shapes = collect_shape_stats(area_crs=area_crs)
    df_build_bus_regions = collect_bus_regions_stats()

    df_only_computational = {
        rname: collect_only_computational(rname)
        for rname in ["build_osm_network", "build_demand_profiles", "build_powerplants"]
    }

    network_dict = {
        network_rule: collect_network_stats(network_rule, scenario_config)
        for network_rule in [
            "base_network",
            "add_electricity",
            "simplify_network",
            "cluster_network",
            "solve_network",
        ]
    }

    # build_renewable_profiles rule
    ren_rule = "build_renewable_profiles"
    renewables_dict = {
        f"{ren_rule}_{tech}": collect_renewable_stats(ren_rule, tech)
        for tech in renewable_config
        if tech in renewable_carriers_config
    }

    # network-related rules
    dict_dfs = {
        "download_osm_data": df_osm_raw,
        "clean_osm_data": df_osm_clean,
        "build_shapes": df_shapes,
        "build_bus_regions": df_build_bus_regions,
        **df_only_computational,
        **renewables_dict,
        **network_dict,
    }

    dict_dfs["total_comp_stats"] = aggregate_computational_stats(
        "total_comp_stats", dict_dfs
    )
    dict_dfs["snakemake_status"] = collect_snakemake_stats(
        "snakemake_status", dict_dfs, renewable_config, renewable_carriers_config
    )

    return dict_dfs


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("make_statistics")

    fp_stats = snakemake.output["stats"]
    scenario = snakemake.params.scenario
    scenario_name = snakemake.config["run"]["name"]

    geo_crs = snakemake.params.crs["geo_crs"]
    metric_crs = snakemake.params.crs["distance_crs"]
    area_crs = snakemake.params.crs["area_crs"]

    renewable = snakemake.params.renewable
    renewable_carriers = snakemake.params.renewable_carriers

    name_index = (
        scenario_name if not scenario_name else "-".join(snakemake.params.countries)
    )

    # create statistics
    stats = calculate_stats(
        scenario,
        renewable,
        renewable_carriers,
        metric_crs=metric_crs,
        area_crs=area_crs,
    )
    stats = pd.concat(stats.values(), axis=1).set_index(pd.Index([name_index]))
    to_csv_nafix(stats, fp_stats)
