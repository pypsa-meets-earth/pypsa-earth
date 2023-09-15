#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 19:01:13 2022

@author: user
"""


import pandas as pd
from prepare_sector_network import get
import os
from pathlib import Path
import country_converter as coco
from helpers import get_conv_factors, aggregate_fuels
# def calc_industry_base(df):

def create_industry_base_totals(df):
    #for country in countries:
    index_mass = df.loc[
        df["Unit"] == "Metric tons,  thousand"
    ].index
    df.loc[index_mass, "Quantity_TWh"] = df.loc[index_mass].apply(
        lambda x: x["Quantity"] * fuels_conv_toTWh[x["Commodity"]], axis=1
    )
    index_energy = df[
        df["Unit"] == "Kilowatt-hours, million"
    ].index
    df.loc[index_energy, "Quantity_TWh"] = df.loc[
        index_energy
    ].apply(lambda x: x["Quantity"] / 1e3, axis=1)

    index_energy_TJ = df[df["Unit"] == "Terajoules"].index
    df.loc[index_energy_TJ, "Quantity_TWh"] = df.loc[
        index_energy_TJ
    ].apply(lambda x: x["Quantity"] / 3600, axis=1)

    index_volume = df[
        df["Unit"] == "Cubic metres, thousand"
    ].index
    df.loc[index_volume, "Quantity_TWh"] = df.loc[
        index_volume
    ].apply(lambda x: x["Quantity"] * fuels_conv_toTWh[x["Commodity"]], axis=1)

    df["carrier"] = df["Commodity"].map(fuel_dict)
    #industry_totals_base = calc_industry_base(df) 

    df_agg = df.groupby(['country', 'carrier', 'Transaction']).agg({'Quantity_TWh': 'sum'}).reset_index()
    industry_totals_base = df_agg.pivot_table(columns='Transaction', index=['country', 'carrier']).fillna(0.0)
    industry_totals_base=industry_totals_base.droplevel(level=0, axis=1)
    return industry_totals_base

if __name__ == "__main__":
    if "snakemake" not in globals():
        from helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        snakemake = mock_snakemake(
            "build_industrial_production_per_country_tomorrow",
            planning_horizons=2030,
            demand="EG",
        )

    config = snakemake.config
    config_ind = snakemake.config["industry"]

    investment_year = int(snakemake.wildcards.planning_horizons)
    demand_sc = snakemake.wildcards.demand

    no_years = int(snakemake.wildcards.planning_horizons) - int(
        snakemake.config["demand_data"]["base_year"]
    )

    all_files = Path("/nfs/home/haz43975/pes_paper/EG/pypsa-earth-sec/data/demand/unsd/data").glob("*.txt") #TODO change path

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
    fuels_conv_toTWh = get_conv_factors("industry")
    # Lists that combine the different fuels in the dataset to the model's carriers
    
    gas_fuels, oil_fuels, biomass_fuels, coal_fuels, heat, electricity = aggregate_fuels("industry")


    fuel_dict = {element: var_name for var_name, element_list in [('gas', gas_fuels), 
                                                                  ('oil', oil_fuels), 
                                                                  ('biomass', biomass_fuels), 
                                                                  ('heat', heat), 
                                                                  ('coal', coal_fuels), 
                                                                  ('electricity', electricity)] for element in element_list}  

    # Fetch country list and demand base year from the config file
    year = snakemake.config["demand_data"]["base_year"]
    countries = snakemake.config["countries"]
    countries = ["EG", "MA"]

    # Filter for the year and country
    df_yr = df[df.Year == year]
    df_yr = df_yr[df_yr.country.isin(countries)]
    df_yr = df_yr[df_yr["Commodity - Transaction"].str.contains("industry")]

    industry_totals_base = create_industry_base_totals(df_yr)

    
    industry_totals_base.to_csv("../data/industry_totals_base.csv")



