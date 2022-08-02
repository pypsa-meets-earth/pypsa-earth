#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 19:01:13 2022

@author: user
"""


import pandas as pd
from prepare_sector_network import get


def industry_prod_tomorrow(production):

    production_tomorrow = production.copy()
    # Steel production

    st_keys = [
        "Integrated steelworks",
        "Electric arc",
    ]  # Columns representing steel prod.
    total_steel = production_tomorrow[st_keys].sum(axis=1)

    if total_steel.sum():
        fraction_prim_st = get(
            config_ind["St_primary_fraction"], investment_year
        )  # Fraction of primary prod. route (integrated)
        fraction_dri = get(
            config_ind["DRI_fraction"], investment_year
        )  # Fraction of the primary to use DRI
        primary_steel = production_tomorrow["Integrated steelworks"].sum()
        fraction_persistent_primary = (
            fraction_prim_st * total_steel.sum() / primary_steel
        )  # Fraction that will remain as primary route (integrated)

        dri = (
            fraction_dri
            * fraction_persistent_primary
            * production_tomorrow["Integrated steelworks"]
        )
        production_tomorrow.insert(2, "DRI + Electric arc", dri)

        not_dri = 1 - fraction_dri
        production_tomorrow["Integrated steelworks"] = (
            not_dri
            * fraction_persistent_primary
            * production_tomorrow["Integrated steelworks"]
        )
        production_tomorrow["Electric arc"] = (
            total_steel
            - production_tomorrow["DRI + Electric arc"]
            - production_tomorrow["Integrated steelworks"]
        )

    # Alimunum production
    keys = ["Aluminium - primary production", "Aluminium - secondary production"]
    total_aluminium = production_tomorrow[keys].sum(axis=1)

    if total_aluminium.sum():
        key_pri = "Aluminium - primary production"
        key_sec = "Aluminium - secondary production"

        al_primary_fraction = get(config_ind["Al_primary_fraction"], investment_year)

        fraction_persistent_primary = (
            al_primary_fraction
            * total_aluminium.sum()
            / production_tomorrow[key_pri].sum()
        )

        production_tomorrow[key_pri] = (
            fraction_persistent_primary * production_tomorrow[key_pri]
        )
        production_tomorrow[key_sec] = total_aluminium - production_tomorrow[key_pri]

    # High Value chemicals
    production_tomorrow["HVC (mechanical recycling)"] = (
        get(config_ind["HVC_mechanical_recycling_fraction"], investment_year)
        * production_tomorrow["HVC"]
    )
    production_tomorrow["HVC (chemical recycling)"] = (
        get(config_ind["HVC_chemical_recycling_fraction"], investment_year)
        * production_tomorrow["HVC"]
    )

    production_tomorrow["HVC"] *= get(
        config_ind["HVC_primary_fraction"], investment_year
    )

    return production_tomorrow


if __name__ == "__main__":
    if "snakemake" not in globals():
        from helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industrial_production_per_country_tomorrow", planning_horizons=2030
        )

    config = snakemake.config
    config_ind = snakemake.config["industry"]

    investment_year = int(snakemake.wildcards.planning_horizons)

    production = pd.read_csv(
        snakemake.input.industrial_production_per_country, index_col=0
    ).filter(config["countries"], axis=0)

    production.at["MA", "Integrated steelworks"] = 1e5
    prod_tomorrow = industry_prod_tomorrow(production)

    prod_tomorrow.to_csv(
        snakemake.output.industrial_production_per_country_tomorrow, float_format="%.2f"
    )
