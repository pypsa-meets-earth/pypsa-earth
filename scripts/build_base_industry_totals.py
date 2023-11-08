#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 19:01:13 2022

@author: user
"""


import os
from pathlib import Path

import country_converter as coco
import pandas as pd
from helpers import aggregate_fuels, get_conv_factors
from prepare_sector_network import get

# def calc_industry_base(df):


def calculate_end_values(df):
    return (1 + df) ** no_years


def create_industry_base_totals(df):
    # Converting values of mass (ktons) to energy (TWh)
    index_mass = df.loc[df["Unit"] == "Metric tons,  thousand"].index
    df.loc[index_mass, "Quantity_TWh"] = df.loc[index_mass].apply(
        lambda x: x["Quantity"] * fuels_conv_toTWh[x["Commodity"]], axis=1
    )

    # Converting values of energy (GWh) to energy (TWh)
    index_energy = df[df["Unit"] == "Kilowatt-hours, million"].index
    df.loc[index_energy, "Quantity_TWh"] = df.loc[index_energy].apply(
        lambda x: x["Quantity"] / 1e3, axis=1
    )

    # Converting values of energy (TJ) to energy (TWh)
    index_energy_TJ = df[df["Unit"] == "Terajoules"].index
    df.loc[index_energy_TJ, "Quantity_TWh"] = df.loc[index_energy_TJ].apply(
        lambda x: x["Quantity"] / 3600, axis=1
    )

    # Converting values of volume (thousand m3) to energy (TWh)
    index_volume = df[df["Unit"] == "Cubic metres, thousand"].index
    df.loc[index_volume, "Quantity_TWh"] = df.loc[index_volume].apply(
        lambda x: x["Quantity"] * fuels_conv_toTWh[x["Commodity"]], axis=1
    )

    df["carrier"] = df["Commodity"].map(fuel_dict)

    # Aggregating and grouping the dataframe
    df_agg = (
        df.groupby(["country", "carrier", "Transaction"])
        .agg({"Quantity_TWh": "sum"})
        .reset_index()
    )
    industry_totals_base = df_agg.pivot_table(
        columns="Transaction", index=["country", "carrier"]
    ).fillna(0.0)
    industry_totals_base = industry_totals_base.droplevel(level=0, axis=1)
    industry_totals_base["other"] = 0

    if include_other:
        # Loop through the columns in the list and sum them if they exist
        for col_name in other_list:
            if col_name in industry_totals_base.columns:
                industry_totals_base["other"] += industry_totals_base[col_name]
                industry_totals_base.drop(col_name, axis=1, inplace=True)
    else:
        industry_totals_base.drop(
            columns=[col for col in other_list if col in industry_totals_base.columns],
            inplace=True,
        )

    industry_totals_base = industry_totals_base.rename(
        columns={"paper, pulp and print": "paper pulp and print"}
    )

    missing_columns = [
        col for col in clean_industry_list if col not in industry_totals_base.columns
    ]

    # Add missing columns with all values set to 0
    for col in missing_columns:
        industry_totals_base[col] = 0

    return industry_totals_base * 1e6  # change from TWh to MWh


if __name__ == "__main__":
    if "snakemake" not in globals():
        from helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        snakemake = mock_snakemake(
            "build_base_industry_totals",
            planning_horizons=2030,
            demand="EG",
        )

    # Loading config file and wild cards

    year = snakemake.config["demand_data"]["base_year"]
    countries = snakemake.config["countries"]
    # countries = ["DE", "US", "EG", "MA", "UA", "UK"]
    # countries = ["EG", "BH"]

    investment_year = int(snakemake.wildcards.planning_horizons)
    demand_sc = snakemake.wildcards.demand
    no_years = int(snakemake.wildcards.planning_horizons) - int(
        snakemake.config["demand_data"]["base_year"]
    )
    include_other = snakemake.config["demand_data"]["other_industries"]

    industry_list = [
        "iron and steel",
        "chemical and petrochemical",
        "non-ferrous metals",
        "non-metallic minerals",
        "transport equipment",
        "machinery",
        "mining and quarrying",
        "food and tobacco",
        "paper, pulp and printing",
        "paper, pulp and print",
        "wood and wood products",
        "textile and leather",
        "construction",
        "not elsewhere specified (industry)",
        "other manuf., const. and non-fuel ind.",
        "other manuf., const. and non-fuel min. ind.",
    ]
    other_list = [
        "other manuf., const. and non-fuel ind.",
        "other manuf., const. and non-fuel min. ind.",
    ]

    clean_industry_list = [
        "iron and steel",
        "chemical and petrochemical",
        "non-ferrous metals",
        "non-metallic minerals",
        "transport equipment",
        "machinery",
        "mining and quarrying",
        "food and tobacco",
        "paper pulp and print",
        "wood and wood products",
        "textile and leather",
        "construction",
        "other",
    ]

    # Loading all energy balance files
    all_files = Path("data/demand/unsd/data").glob("*.txt")  # TODO change path

    # Create a dataframe from all downloaded files
    df = pd.concat(
        (pd.read_csv(f, encoding="utf8", sep=";") for f in all_files), ignore_index=True
    )

    # Split 'Commodity', 'Transaction' column to two
    df[["Commodity", "Transaction", "extra"]] = df["Commodity - Transaction"].str.split(
        " - ", expand=True
    )

    # Remove fill na in Transaction column
    df["Transaction"] = df["Transaction"].fillna("NA")
    df["Transaction"] = df["Transaction"].str.lower()
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
    fuels_conv_toTWh = get_conv_factors("industry")

    # Lists that combine the different fuels in the dataset to the model's carriers

    # Fetch the fuel categories from the helpers script
    (
        gas_fuels,
        oil_fuels,
        biomass_fuels,
        coal_fuels,
        heat,
        electricity,
    ) = aggregate_fuels("industry")

    # Create fuel dictionary to use for mapping all fuels to the pypsa representative fuels
    fuel_dict = {
        element: var_name
        for var_name, element_list in [
            ("gas", gas_fuels),
            ("oil", oil_fuels),
            ("biomass", biomass_fuels),
            ("heat", heat),
            ("coal", coal_fuels),
            ("electricity", electricity),
        ]
        for element in element_list
    }

    # Filter for the year and country
    df_yr = df[df.Year == year]
    df_yr = df_yr[df_yr.country.isin(countries)]
    # df_yr = df_yr[df_yr["Commodity - Transaction"].str.contains("industry")]
    df_yr = df_yr[
        df_yr["Transaction"].str.lower().str.contains("|".join(industry_list))
    ]

    df_yr = df_yr[
        ~df_yr.Transaction.str.lower().str.contains(
            "consumption by manufacturing, construction and non-fuel"
        )
    ]
    df_yr["Transaction"] = df_yr["Transaction"].str.replace(" industry", "")

    df_yr["Transaction"] = df_yr["Transaction"].str.rstrip()
    df_yr["Transaction"] = df_yr["Transaction"].str.replace("consumption by ", "")

    # df_yr[df_yr['Commodity - Transaction'].str.lower().str.contains('|'.join(industry_list))]
    # df_yr "consumption by manufacturing, construction and non-fuel"
    # Create the industry totals file
    industry_totals_base = create_industry_base_totals(df_yr)

    # Export the industry totals dataframe
    industry_totals_base.to_csv(snakemake.output["base_industry_totals"])

    # print("end")
