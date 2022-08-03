#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 21:18:06 2022

@author: user
"""

from itertools import product

import numpy as np
import pandas as pd

# # map JRC/our sectors to hotmaps sector, where mapping exist
# sector_mapping = {
#     'Electric arc': 'Iron and steel',
#     'Integrated steelworks': 'Iron and steel',
#     'DRI + Electric arc': 'Iron and steel',
#     'Ammonia': 'Chemical industry',
#     'HVC': 'Chemical industry',
#     'HVC (mechanical recycling)': 'Chemical industry',
#     'HVC (chemical recycling)': 'Chemical industry',
#     'Methanol': 'Chemical industry',
#     'Chlorine': 'Chemical industry',
#     'Other chemicals': 'Chemical industry',
#     'Pharmaceutical products etc.': 'Chemical industry',
#     'Cement': 'Cement',
#     'Ceramics & other NMM': 'Non-metallic mineral products',
#     'Glass production': 'Glass',
#     'Pulp production': 'Paper and printing',
#     'Paper production': 'Paper and printing',
#     'Printing and media reproduction': 'Paper and printing',
#     'Alumina production': 'Non-ferrous metals',
#     'Aluminium - primary production': 'Non-ferrous metals',
#     'Aluminium - secondary production': 'Non-ferrous metals',
#     'Other non-ferrous metals': 'Non-ferrous metals',
# }


def country_to_nodal(industrial_production, keys):

    keys["country"] = keys.index.str[:2]  # TODO 2digit_3_digit adaptation needed

    nodal_production = pd.DataFrame(
        index=keys.index, columns=industrial_production.columns, dtype=float
    )

    countries = keys.country.unique()
    sectors = industrial_production.columns

    for country, sector in product(countries, sectors):

        buses = keys.index[keys.country == country]

        if sector not in dist_keys.columns or dist_keys[sector].sum() == 0:
            mapping = "population"
        else:
            mapping = sector

        key = keys.loc[buses, mapping]
        nodal_production.loc[buses, sector] = (
            industrial_production.at[country, sector] * key
        )

    return nodal_production


if __name__ == "__main__":
    if "snakemake" not in globals():
        from helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industry_demand", simpl="", clusters=9077, planning_horizons="2030"
        )

    # Load production per country tomorrow
    prod_tom_path = snakemake.input.industrial_production_per_country_tomorrow
    production_tom = pd.read_csv(prod_tom_path, header=0, index_col=0)

    # Load distribution keys
    keys_path = snakemake.input.industrial_distribution_key
    dist_keys = pd.read_csv(keys_path, index_col=0)

    # material demand per node and industry (kton/a)
    nodal_production = country_to_nodal(production_tom, dist_keys)

    # Transfromation key to map the material demand to the corresponding carrier demand
    industry_sector_ratios = pd.read_csv(
        snakemake.input.industry_sector_ratios, index_col=0
    )

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

    # energy demand today to get current electricity #TODO
    prod_tod_path = snakemake.input.industrial_production_per_country
    production_tod = pd.read_csv(prod_tod_path, header=0, index_col=0).filter(
        snakemake.config["countries"], axis=0
    )

    nodal_production_tod = country_to_nodal(production_tod, dist_keys)
    nodal_production_tod = nodal_production_tod.reindex(
        columns=list(industry_sector_ratios.columns), fill_value=0
    )
    nodal_energy_today = nodal_production_tod.dot(industry_sector_ratios.T)
    nodal_energy_today.rename(columns=rename_sectors, inplace=True)
    nodal_energy_today *= 0.001

    nodal_df["current electricity"] = nodal_energy_today["electricity"]

    nodal_df.index.name = "TWh/a (MtCO2/a)"

    nodal_df.to_csv(
        snakemake.output.industrial_energy_demand_per_node, float_format="%.2f"
    )

    # nodal_demand = build_nodal_industrial_energy_demand(demand_ct_td, dist_keys)
    # nodal_demand.to_csv(snakemake.output.industrial_energy_demand_per_node_today)
