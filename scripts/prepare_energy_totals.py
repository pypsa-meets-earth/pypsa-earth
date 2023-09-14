# -*- coding: utf-8 -*-
import glob
import os
import sys
from io import BytesIO
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile

import country_converter as coco
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import py7zr
import requests
from helpers import sets_path_to_root, three_2_two_digits_country


def get(item, investment_year=None):
    """Check whether item depends on investment year"""
    if isinstance(item, dict):
        return item[investment_year]
    else:
        return item


def calculate_end_values(df):
    return (1 + df) ** no_years


if __name__ == "__main__":
    if "snakemake" not in globals():
        from helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        snakemake = mock_snakemake(
            "prepare_energy_totals",
            simpl="",
            clusters=32,
            demand="EG",
            planning_horizons=2030,
        )
        sets_path_to_root("pypsa-earth-sec")

    countries = snakemake.config["countries"]
    # countries = ["NG", "BJ"]
    investment_year = int(snakemake.wildcards.planning_horizons)
    demand_sc = snakemake.wildcards.demand  # loading the demand scenrario wildcard

    base_energy_totals = pd.read_csv("data/energy_totals_base.csv", index_col=0)
    growth_factors_cagr = pd.read_csv(
        "data/demand/growth_factors_cagr.csv", index_col=0
    )
    efficiency_gains_cagr = pd.read_csv(
        "data/demand/efficiency_gains_cagr.csv", index_col=0
    )
    fuel_shares = pd.read_csv("data/demand/fuel_shares.csv", index_col=0)
    district_heating = pd.read_csv("data/demand/district_heating.csv", index_col=0)

    no_years = int(snakemake.wildcards.planning_horizons) - int(
        snakemake.config["demand_data"]["base_year"]
    )
    growth_factors = calculate_end_values(growth_factors_cagr)
    efficiency_gains = calculate_end_values(efficiency_gains_cagr)
    # efficiency_gains = efficiency_gains[efficiency_gains.index.isin(countries)]
    # fuel_shares = fuel_shares[fuel_shares.index.isin(countries)]
    # district_heating = district_heating[district_heating.index.isin(countries)]

    options = snakemake.config["sector"]

    fuel_cell_share = get(
        options["land_transport_fuel_cell_share"],
        demand_sc + "_" + str(investment_year),
    )
    electric_share = get(
        options["land_transport_electric_share"], demand_sc + "_" + str(investment_year)
    )

    hydrogen_shipping_share = get(
        options["shipping_hydrogen_share"], demand_sc + "_" + str(investment_year)
    )

    energy_totals = (
        base_energy_totals
        * efficiency_gains.loc[countries]
        * growth_factors.loc[countries]
    )

    # Residential
    efficiency_heat_oil_to_elec = snakemake.config["sector"][
        "efficiency_heat_oil_to_elec"
    ]
    efficiency_heat_biomass_to_elec = snakemake.config["sector"][
        "efficiency_heat_biomass_to_elec"
    ]
    efficiency_heat_gas_to_elec = snakemake.config["sector"][
        "efficiency_heat_gas_to_elec"
    ]

    energy_totals["electricity residential space"] = (
        base_energy_totals["total residential space"]
        + (
            fuel_shares["biomass to elec heat share"]
            * fuel_shares["biomass residential heat share"]
            * (fuel_shares["space to water heat share"])
            * base_energy_totals["residential biomass"]
            * efficiency_heat_biomass_to_elec
        )
        + (
            fuel_shares["oil to elec heat share"]
            * fuel_shares["oil residential heat share"]
            * (fuel_shares["space to water heat share"])
            * base_energy_totals["residential oil"]
            * efficiency_heat_oil_to_elec
        )
        + (
            fuel_shares["gas to elec heat share"]
            * fuel_shares["gas residential heat share"]
            * (fuel_shares["space to water heat share"])
            * base_energy_totals["residential gas"]
            * efficiency_heat_gas_to_elec
        )
    )

    energy_totals["electricity residential water"] = (
        base_energy_totals["total residential water"]
        + (
            fuel_shares["biomass to elec heat share"]
            * fuel_shares["biomass residential heat share"]
            * (1 - fuel_shares["space to water heat share"])
            * base_energy_totals["residential biomass"]
            * efficiency_heat_biomass_to_elec
        )
        + (
            fuel_shares["oil to elec heat share"]
            * fuel_shares["oil residential heat share"]
            * (1 - fuel_shares["space to water heat share"])
            * base_energy_totals["residential oil"]
            * efficiency_heat_oil_to_elec
        )
        + (
            fuel_shares["gas to elec heat share"]
            * fuel_shares["gas residential heat share"]
            * (1 - fuel_shares["space to water heat share"])
            * base_energy_totals["residential gas"]
            * efficiency_heat_gas_to_elec
        )
    )

    energy_totals["residential heat oil"] = (
        base_energy_totals["residential oil"]
        * fuel_shares["oil residential heat share"]
        * (1 - fuel_shares["oil to elec heat share"])
        * efficiency_gains["residential heat oil"]
        * growth_factors["residential heat oil"]
    )

    energy_totals["residential oil"] = (
        base_energy_totals["residential oil"]
        * (1 - fuel_shares["oil residential heat share"])
        * (1 - fuel_shares["oil to elec share"])
        * efficiency_gains["residential oil"]
        * growth_factors["residential oil"]
    )

    energy_totals["residential heat biomass"] = (
        base_energy_totals["residential biomass"]
        * fuel_shares["biomass residential heat share"]
        * (1 - fuel_shares["biomass to elec heat share"])
        * efficiency_gains["residential heat biomass"]
        * growth_factors["residential heat biomass"]
    )

    energy_totals["residential biomass"] = (
        base_energy_totals["residential biomass"]
        * (1 - fuel_shares["biomass residential heat share"])
        * (1 - fuel_shares["biomass to elec share"])
        * efficiency_gains["residential biomass"]
        * growth_factors["residential biomass"]
    )

    energy_totals["residential heat gas"] = (
        base_energy_totals["residential gas"]
        * fuel_shares["gas residential heat share"]
        * (1 - fuel_shares["gas to elec heat share"])
        * efficiency_gains["residential heat gas"]
        * growth_factors["residential heat gas"]
    )

    energy_totals["residential gas"] = (
        base_energy_totals["residential gas"]
        * (1 - fuel_shares["gas residential heat share"])
        * (1 - fuel_shares["gas to elec share"])
        * efficiency_gains["residential gas"]
        * growth_factors["residential gas"]
    )

    energy_totals["total residential space"] = energy_totals[
        "electricity residential space"
    ] + (
        energy_totals["residential heat oil"]
        + energy_totals["residential heat biomass"]
        + energy_totals["residential heat gas"]
    ) * (
        fuel_shares["space to water heat share"]
    )

    energy_totals["total residential water"] = energy_totals[
        "electricity residential water"
    ] + (
        energy_totals["residential heat oil"]
        + energy_totals["residential heat biomass"]
        + energy_totals["residential heat gas"]
    ) * (
        1 - fuel_shares["space to water heat share"]
    )

    energy_totals["electricity residential"] = (
        energy_totals["electricity residential"]
        + (
            fuel_shares["oil to elec share"]
            * (1 - fuel_shares["oil residential heat share"])
            * base_energy_totals["residential oil"]
        )
        + (
            fuel_shares["biomass to elec share"]
            * (1 - fuel_shares["biomass residential heat share"])
            * base_energy_totals["residential biomass"]
        )
        + (
            fuel_shares["gas to elec share"]
            * (1 - fuel_shares["gas residential heat share"])
            * base_energy_totals["residential gas"]
        )
    )

    # Road
    energy_totals["total road"] = (
        (1 - fuel_cell_share - electric_share)
        * efficiency_gains["total road ice"]
        * base_energy_totals["total road"]
        + fuel_cell_share
        * efficiency_gains["total road fcev"]
        * base_energy_totals["total road"]
        + electric_share
        * efficiency_gains["total road ev"]
        * base_energy_totals["total road"]
    ) * growth_factors["total road"]

    # Navigation
    energy_totals["total domestic navigation"] = (
        (1 - hydrogen_shipping_share)
        * efficiency_gains["total navigation oil"]
        * base_energy_totals["total domestic navigation"]
        + hydrogen_shipping_share
        * efficiency_gains["total navigation hydrogen"]
        * base_energy_totals["total domestic navigation"]
    ) * growth_factors["total domestic navigation"]

    energy_totals["total international navigation"] = (
        (1 - hydrogen_shipping_share)
        * efficiency_gains["total navigation oil"]
        * base_energy_totals["total international navigation"]
        + hydrogen_shipping_share
        * efficiency_gains["total navigation hydrogen"]
        * base_energy_totals["total international navigation"]
    ) * growth_factors["total international navigation"]

    energy_totals["district heat share"] = district_heating["current"]

    energy_totals["electricity services space"] = 0
    energy_totals["electricity services water"] = 0

    energy_totals = energy_totals.dropna(axis=1, how="all")

    energy_totals.fillna(0).to_csv(snakemake.output.energy_totals)
