# -*- coding: utf-8 -*-
import glob
import logging
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

_logger = logging.getLogger(__name__)

pd.options.mode.chained_assignment = None


def calc_sector(sector):
    for country in countries:
        # print(country, sector)
        df_co = df_yr[df_yr.country == country]

        if sector != "navigation":
            df_sector = df_co.loc[
                df["Commodity - Transaction"].str.lower().str.contains(sector)
            ]
            # assert df_yr[df_yr["Commodity - Transaction"].str.contains(sector)]["Unit"].unique() == 'Metric tons,  thousand', "Not all quantities have the expected unit: {}".format(expected_unit)
        else:
            df_sector = df_co.loc[
                (df["Commodity - Transaction"].str.lower().str.contains(sector))
                | (
                    df["Commodity - Transaction"]
                    .str.lower()
                    .str.contains("marine bunkers")
                )
            ]

        if df_sector.empty:
            if sector == "consumption by households":
                energy_totals_base.at[country, "electricity residential"] = np.NaN
                energy_totals_base.at[country, "residential oil"] = np.NaN
                energy_totals_base.at[country, "residential biomass"] = np.NaN
                energy_totals_base.at[country, "residential gas"] = np.NaN
                energy_totals_base.at[country, "total residential space"] = np.NaN
                energy_totals_base.at[country, "total residential water"] = np.NaN

            elif sector == "services":
                energy_totals_base.at[country, "services electricity"] = np.NaN
                energy_totals_base.at[country, "services oil"] = np.NaN
                energy_totals_base.at[country, "services biomass"] = np.NaN
                energy_totals_base.at[country, "services gas"] = np.NaN
                energy_totals_base.at[country, "total services space"] = np.NaN
                energy_totals_base.at[country, "total services water"] = np.NaN

            elif sector == "road":
                energy_totals_base.at[country, "total road"] = np.NaN

            elif sector == "agriculture":
                energy_totals_base.at[country, "agriculture electricity"] = np.NaN
                energy_totals_base.at[country, "agriculture oil"] = np.NaN
                energy_totals_base.at[country, "agriculture biomass"] = np.NaN
                # energy_totals_base.at[country, "electricity rail"] = np.NaN

            elif sector == "rail":
                energy_totals_base.at[country, "total rail"] = np.NaN
                energy_totals_base.at[country, "electricity rail"] = np.NaN

            elif sector == "aviation":
                energy_totals_base.at[country, "total international aviation"] = np.NaN
                energy_totals_base.at[country, "total domestic aviation"] = np.NaN

            elif sector == "navigation":
                energy_totals_base.at[
                    country, "total international navigation"
                ] = np.NaN
                energy_totals_base.at[country, "total domestic navigation"] = np.NaN

            _logger.warning("No data for " + country + " in the sector " + sector + ".")

        else:
            index_mass = df_sector.loc[
                df_sector["Unit"] == "Metric tons,  thousand"
            ].index
            df_sector.loc[index_mass, "Quantity_TWh"] = df_sector.loc[index_mass].apply(
                lambda x: x["Quantity"] * fuels_conv_toTWh[x["Commodity"]], axis=1
            )

            index_energy = df_sector[
                df_sector["Unit"] == "Kilowatt-hours, million"
            ].index
            df_sector.loc[index_energy, "Quantity_TWh"] = df_sector.loc[
                index_energy
            ].apply(lambda x: x["Quantity"] / 1e3, axis=1)

            index_energy_TJ = df_sector[df_sector["Unit"] == "Terajoules"].index
            df_sector.loc[index_energy_TJ, "Quantity_TWh"] = df_sector.loc[
                index_energy_TJ
            ].apply(lambda x: x["Quantity"] / 3600, axis=1)

            index_volume = df_sector[
                df_sector["Unit"] == "Cubic metres, thousand"
            ].index
            df_sector.loc[index_volume, "Quantity_TWh"] = df_sector.loc[
                index_volume
            ].apply(lambda x: x["Quantity"] * fuels_conv_toTWh[x["Commodity"]], axis=1)

            sectors_dfs[sector] = df_sector.copy()

            if sector == "consumption by households":
                energy_totals_base.at[country, "electricity residential"] = round(
                    df_sector[
                        (df_sector.Commodity == "Electricity")
                        | df_sector.Commodity.isin(other_fuels)
                    ].Quantity_TWh.sum(),
                    4,
                )
                energy_totals_base.at[country, "residential oil"] = round(
                    df_sector[df_sector.Commodity.isin(oil_fuels)].Quantity_TWh.sum(), 4
                )
                energy_totals_base.at[country, "residential biomass"] = round(
                    df_sector[
                        df_sector.Commodity.isin(biomass_fuels)
                    ].Quantity_TWh.sum(),
                    4,
                )
                energy_totals_base.at[country, "residential gas"] = round(
                    df_sector[df_sector.Commodity.isin(gas_fuels)].Quantity_TWh.sum(), 4
                )
                energy_totals_base.at[country, "total residential space"] = (
                    round(
                        df_sector[df_sector.Commodity.isin(heat)].Quantity_TWh.sum(), 4
                    )
                    * snakemake.config["sector"]["space_heat_share"]
                )
                energy_totals_base.at[country, "total residential water"] = round(
                    df_sector[df_sector.Commodity.isin(heat)].Quantity_TWh.sum(), 4
                ) * (1 - snakemake.config["sector"]["space_heat_share"])

            elif sector == "services":
                energy_totals_base.at[country, "services electricity"] = round(
                    df_sector[
                        (df_sector.Commodity == "Electricity")
                        | df_sector.Commodity.isin(other_fuels)
                    ].Quantity_TWh.sum(),
                    4,
                )
                energy_totals_base.at[country, "services oil"] = round(
                    df_sector[df_sector.Commodity.isin(oil_fuels)].Quantity_TWh.sum(), 4
                )
                energy_totals_base.at[country, "services biomass"] = round(
                    df_sector[
                        df_sector.Commodity.isin(biomass_fuels)
                    ].Quantity_TWh.sum(),
                    4,
                )
                energy_totals_base.at[country, "services gas"] = round(
                    df_sector[df_sector.Commodity.isin(gas_fuels)].Quantity_TWh.sum(), 4
                )
                energy_totals_base.at[country, "total services space"] = (
                    round(
                        df_sector[df_sector.Commodity.isin(heat)].Quantity_TWh.sum(), 4
                    )
                    * snakemake.config["sector"]["space_heat_share"]
                )
                energy_totals_base.at[country, "total services water"] = round(
                    df_sector[df_sector.Commodity.isin(heat)].Quantity_TWh.sum(), 4
                ) * (1 - snakemake.config["sector"]["space_heat_share"])

            elif sector == "road":
                energy_totals_base.at[country, "total road"] = round(
                    df_sector.Quantity_TWh.sum(), 4
                )

            elif sector == "agriculture":
                energy_totals_base.at[country, "agriculture electricity"] = round(
                    df_sector[
                        (df_sector.Commodity == "Electricity")
                        | df_sector.Commodity.isin(other_fuels)
                    ].Quantity_TWh.sum(),
                    4,
                )
                energy_totals_base.at[country, "agriculture oil"] = round(
                    df_sector[df_sector.Commodity.isin(oil_fuels)].Quantity_TWh.sum(), 4
                )
                energy_totals_base.at[country, "agriculture biomass"] = round(
                    df_sector[
                        df_sector.Commodity.isin(biomass_fuels)
                    ].Quantity_TWh.sum(),
                    4,
                )
                # energy_totals_base.at[country, "electricity rail"] = round(df_house[(df_house.Commodity=="Electricity")].Quantity_TWh.sum(), 4)

            elif sector == "rail":
                energy_totals_base.at[country, "total rail"] = round(
                    df_sector[
                        (df_sector.Commodity == "Gas Oil/ Diesel Oil")
                        | (df_sector.Commodity == "Biodiesel")
                        | (df_sector.Commodity == "Electricity")
                    ].Quantity_TWh.sum(),
                    4,
                )
                energy_totals_base.at[country, "electricity rail"] = round(
                    df_sector[
                        (df_sector.Commodity == "Electricity")
                    ].Quantity_TWh.sum(),
                    4,
                )

            elif sector == "aviation":
                energy_totals_base.at[country, "total international aviation"] = round(
                    df_sector[
                        (df_sector.Commodity == "Kerosene-type Jet Fuel")
                        & (df_sector.Transaction == "International aviation bunkers")
                    ].Quantity_TWh.sum(),
                    4,
                )
                energy_totals_base.at[country, "total domestic aviation"] = round(
                    df_sector[
                        (df_sector.Commodity == "Kerosene-type Jet Fuel")
                        & (df_sector.Transaction == "Consumption by domestic aviation")
                    ].Quantity_TWh.sum(),
                    4,
                )

            elif sector == "navigation":
                energy_totals_base.at[
                    country, "total international navigation"
                ] = round(
                    df_sector[
                        df_sector.Transaction == "International marine bunkers"
                    ].Quantity_TWh.sum(),
                    4,
                )
                energy_totals_base.at[country, "total domestic navigation"] = round(
                    df_sector[
                        df_sector.Transaction == "Consumption by domestic navigation"
                    ].Quantity_TWh.sum(),
                    4,
                )

            else:
                print("wrong sector")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        snakemake = mock_snakemake(
            "build_base_energy_totals",
            simpl="",
            clusters=10,
            demand="DF",
            planning_horizons=2030,
        )
        sets_path_to_root("pypsa-earth-sec")

    energy_stat_database = pd.read_excel(
        snakemake.input.unsd_paths, index_col=0, header=0
    )  # pd.read_excel("/nfs/home/haz43975/pypsa-earth-sec/scripts/Energy_Statistics_Database.xlsx"

    # Load the links and make a dictionary
    df = energy_stat_database.copy()
    df = df.dropna(axis=0, subset=["Link"])
    df = df.to_dict("dict")
    d = df["Link"]

    if snakemake.config["demand_data"]["update_data"]:
        # Delete and existing files to avoid duplication and double counting

        files = glob.glob("data/demand/unsd/data/*.txt")
        for f in files:
            os.remove(f)

        # Feed the dictionary of links to the for loop, download and unzip all files
        for key, value in d.items():
            zipurl = value

            with urlopen(zipurl) as zipresp:
                with ZipFile(BytesIO(zipresp.read())) as zfile:
                    zfile.extractall("data/demand/unsd/data")

                    path = "data/demand/unsd/data"

    # Get the files from the path provided in the OP
    all_files = list(Path("data/demand/unsd/data").glob("*.txt"))

    # Create a dataframe from all downloaded files
    df = pd.concat(
        (pd.read_csv(f, encoding="utf8", sep=";") for f in all_files), ignore_index=True
    )

    # Split 'Commodity', 'Transaction' column to two
    df[["Commodity", "Transaction", "extra"]] = df["Commodity - Transaction"].str.split(
        " - ", expand=True
    )

    # Remove Foootnote and Estimate from 'Commodity - Transaction' column
    df = df.loc[df["Commodity - Transaction"] != "Footnote"]
    df = df.loc[df["Commodity - Transaction"] != "Estimate"]

    # Create a column with iso2 country code
    cc = coco.CountryConverter()
    Country = pd.Series(df["Country or Area"])

    df["country"] = cc.pandas_convert(series=Country, to="ISO2", not_found="not found")

    # remove countries or areas that have no iso2 such as former countries names
    df = df.loc[df["country"] != "not found"]

    # Convert country column that contains lists for some country names that are identified with more than one country.
    df["country"] = df["country"].astype(str)

    # Remove all iso2 conversions for some country names that are identified with more than one country.
    df = df[~df.country.str.contains(",", na=False)].reset_index(drop=True)

    # Create a dictionary with all the conversion factors from ktons or m3 to TWh based on https://unstats.un.org/unsd/energy/balance/2014/05.pdf
    fuels_conv_toTWh = {
        "Gas Oil/ Diesel Oil": 0.01194,
        "Motor Gasoline": 0.01230,
        "Kerosene-type Jet Fuel": 0.01225,
        "Aviation gasoline": 0.01230,
        "Biodiesel": 0.01022,
        "Natural gas liquids": 0.01228,
        "Biogasoline": 0.007444,
        "Bitumen": 0.01117,
        "Fuel oil": 0.01122,
        "Liquefied petroleum gas (LPG)": 0.01313,
        "Liquified Petroleum Gas (LPG)": 0.01313,
        "Lubricants": 0.01117,
        "Naphtha": 0.01236,
        "Fuelwood": 0.00254,
        "Charcoal": 0.00819,
        "Patent fuel": 0.00575,
        "Brown coal briquettes": 0.00575,
        "Hard coal": 0.007167,
        "Other bituminous coal": 0.005556,
        "Anthracite": 0.005,
        "Peat": 0.00271,
        "Peat products": 0.00271,
        "Lignite": 0.003889,
        "Brown coal": 0.003889,
        "Sub-bituminous coal": 0.005555,
        "Coke-oven coke": 0.0002778,
        "Coke oven coke": 0.0002778,
        "Coke Oven Coke": 0.0002778,
        "Gasoline-type jet fuel": 0.01230,
        "Conventional crude oil": 0.01175,
        "Brown Coal Briquettes": 0.00575,
        "Refinery Gas": 0.01375,
        "Petroleum coke": 0.009028,
        "Coking coal": 0.007833,
        "Peat Products": 0.00271,
        "Petroleum Coke": 0.009028,
    }

    # Fetch country list and demand base year from the config file
    year = snakemake.config["demand_data"]["base_year"]
    countries = snakemake.config["countries"]
    # countries = ["NG", "BJ"]

    # Filter for the year and country
    df_yr = df[df.Year == year]
    df_yr = df_yr[df_yr.country.isin(countries)]

    # Create an empty dataframe for energy_totals_base
    energy_totals_cols = pd.read_csv("data/energy_totals_DF_2030.csv").columns
    energy_totals_base = pd.DataFrame(columns=energy_totals_cols, index=countries)

    # Lists that combine the different fuels in the dataset to the model's carriers
    oil_fuels = [
        "Patent fuel",
        "Gas Oil/ Diesel Oil",
        "Motor Gasoline",
        "Liquefied petroleum gas (LPG)",
    ]
    gas_fuels = ["Natural gas (including LNG)", "Gasworks Gas", "Natural gas (including LNG)"]
    biomass_fuels = ["Biodiesel", "Biogases", "Fuelwood"]
    other_fuels = ["Charcoal", "Brown coal briquettes", "Other bituminous coal"]
    heat = ["Heat", "Direct use of geothermal heat", "Direct use of solar thermal heat"]

    # Create a dictionary to save the data if need to be checked
    sectors_dfs = {}

    # Run the function that processes the data for all the sectors
    sectors = [
        "consumption by households",
        "road",
        "rail",
        "aviation",
        "navigation",
        "agriculture",
        "services",
    ]
    for sector in sectors:
        calc_sector(sector)

    # Export the base energy totals file
    energy_totals_base.to_csv(snakemake.output.energy_totals_base)
