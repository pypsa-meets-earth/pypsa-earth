# -*- coding: utf-8 -*-
"""Build industrial energy demand per node."""

import os
from itertools import product

import pandas as pd

sector_mapping = {"Cement": "Cement"}


def build_nodal_industrial_production():

    fn = (
        snakemake.input.industrial_production_per_country_tomorrow
    )  # TODO implement HyPAT data
    industrial_production = pd.read_csv(fn, index_col=0)

    fn = snakemake.input.industrial_distribution_key
    keys = pd.read_csv(fn, index_col=0)
    keys["country"] = keys.index.str[:2]

    nodal_production = pd.DataFrame(
        index=keys.index,
        columns=industrial_production.columns,  # TODO need hypat
        dtype=float,
    )

    countries = keys.country.unique()
    sectors = industrial_production.columns  # TODO need hypat

    for country, sector in product(countries, sectors):

        buses = keys.index[keys.country == country]
        mapping = sector_mapping.get(sector, "population")

        key = keys.loc[buses, mapping]
        nodal_production.loc[buses, sector] = (
            industrial_production.at[country, sector] * key
        )  # TODO need hypat

    # nodal_production.to_csv(snakemake.output.industrial_production_per_node)
    return nodal_production


if __name__ == "__main__":
    if "snakemake" not in globals():
        from helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        snakemake = mock_snakemake(
            "build_industrial_energy_demand_per_node",
            simpl="",
            clusters=48,
            planning_horizons=2030,
        )

    # build_nodal_industrial_production()

    # import EU ratios df as csv
    fn = snakemake.input.industry_sector_ratios
    industry_sector_ratios = pd.read_csv(fn, index_col=0)

    # material demand per node and industry (kton/a)
    # fn = snakemake.input.industrial_production_per_node
    # nodal_production = pd.read_csv(fn, index_col=0)
    fn = build_nodal_industrial_production()
    nodal_production = pd.DataFrame(fn)

    # energy demand today to get current electricity
    fn = snakemake.input.industrial_energy_demand_per_node_today
    nodal_today = pd.read_csv(fn, index_col=0)

    # final energy consumption per node and industry (TWh/a)
    nodal_df = nodal_production.dot(industry_sector_ratios.T)

    # convert GWh to TWh and ktCO2 to MtCO2
    nodal_df *= 0.001

    rename_sectors = {
        "elec": "electricity",
        "biomass": "solid biomass",
        "heat": "low-temperature heat",
    }
    nodal_df.rename(columns=rename_sectors, inplace=True)

    nodal_df["current electricity"] = nodal_today["electricity"]

    nodal_df.index.name = "TWh/a (MtCO2/a)"

    fn = snakemake.output.industrial_energy_demand_per_node
    nodal_df.to_csv(fn, float_format="%.2f")
