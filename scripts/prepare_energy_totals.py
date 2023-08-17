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


if __name__ == "__main__":
    if "snakemake" not in globals():
        from helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        snakemake = mock_snakemake(
            "prepare_energy_totals",
            simpl="",
            clusters=10,
            demand="DF",
            planning_horizons=2030,
        )
        sets_path_to_root("pypsa-earth-sec")

    base_energy_totals = pd.read_csv("data/energy_totals_base.csv", index_col=0)
    growth_factors = pd.read_csv("data/demand/growth_factors.csv", index_col=0)
    efficiency_gains = pd.read_csv("data/demand/efficiency_gains.csv", index_col=0)

    energy_totals = base_energy_totals * efficiency_gains * growth_factors

    options = snakemake.config["sector"]

    investment_year = int(snakemake.wildcards.planning_horizons)
    demand_sc = snakemake.wildcards.demand  # loading the demand scenrario wildcard

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

    energy_totals["total road"] = (
        (1 - fuel_cell_share + electric_share)
        * efficiency_gains["total road ice"]
        * base_energy_totals["total road"]
        + fuel_cell_share
        * efficiency_gains["total road fcev"]
        * base_energy_totals["total road"]
        + electric_share
        * efficiency_gains["total road ev"]
        * base_energy_totals["total road"]
    )

    energy_totals["total domestic navigation"] = (
        1 - hydrogen_shipping_share
    ) * efficiency_gains["total navigation oil"] * base_energy_totals[
        "total domestic navigation"
    ] + hydrogen_shipping_share * efficiency_gains[
        "total navigation hydrogen"
    ] * base_energy_totals[
        "total domestic navigation"
    ]

    energy_totals["total international navigation"] = (
        1 - hydrogen_shipping_share
    ) * efficiency_gains["total navigation oil"] * base_energy_totals[
        "total international navigation"
    ] + hydrogen_shipping_share * efficiency_gains[
        "total navigation hydrogen"
    ] * base_energy_totals[
        "total international navigation"
    ]

    energy_totals["district heat share"] = 0
    energy_totals["electricity residential water"] = 0
    energy_totals["electricity residential space"] = 0
    energy_totals["electricity services space"] = 0
    energy_totals["electricity services water"] = 0

    energy_totals = energy_totals.dropna(axis=1, how="all")

    energy_totals.to_csv(snakemake.output.energy_totals)
