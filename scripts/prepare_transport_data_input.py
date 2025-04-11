# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
import logging
import os
import shutil
from pathlib import Path

import country_converter as coco
import numpy as np
import pandas as pd
from _helpers import BASE_DIR

logger = logging.getLogger(__name__)

def add_iso2_country_code(df):
    """
    Converts 'Country' names to ISO2 codes in a new 'country' column.
    Cleans DataFrame by removing rows with invalid 'country' values.
    """

    cc = coco.CountryConverter()
    df.loc[:, "country"] = cc.pandas_convert(
        series=pd.Series(df["Country"]), to="ISO2", not_found="not found"
    )

    df = df[df.country != "not found"]

    # Drop region names where country column contains list of countries
    df = df.loc[df.country.apply(lambda x: isinstance(x, str)),:]

    return df

def download_number_of_vehicles():
    """
    Downloads and returns the number of registered vehicles as tabular data 
    from the Global Health Observatory (GHO) repository data and from Wikipedia.

    The csv data from the WHO website is imported
    from 'https://apps.who.int/gho/data/node.main.A995'.
    A few countries are missing in the WHO list (e.g. South Africa, Algeria).
    Therefore, the number of vehicles per country table from Wikipedia
    is also imported for completion (prio 2):
    'https://en.wikipedia.org/wiki/List_of_countries_and_territories_by_motor_vehicles_per_capita'.
    """
    
    def _clean_data(df):
        df = df.dropna(subset=["number cars"])
        df.loc[:,"number cars"] = df.loc[:,"number cars"].astype(int)
        return df#[["Country", "number cars"]]
    
    storage_options = {"User-Agent": "Mozilla/5.0"}
    url = "https://apps.who.int/gho/athena/data/GHO/RS_194?filter=COUNTRY:*&ead=&x-sideaxis=COUNTRY;YEAR;DATASOURCE&x-topaxis=GHO&profile=crosstable&format=csv"
    try:
        vehicles_gho = pd.read_csv(
            url, storage_options=storage_options, encoding="utf8"
        )
        print("File read successfully.")
    except Exception as e:
        logger.warning("Failed to read the file. Falling back on hard-coded data:",e)
        return pd.DataFrame()

    vehicles_gho = vehicles_gho.rename(
        columns={
            "Countries, territories and areas": "Country",
            "Number of registered vehicles": "number cars",
        }
    )

    vehicles_gho = add_iso2_country_code(vehicles_gho)

    vehicles_gho["number cars"] = (
        vehicles_gho["number cars"].str.replace(" ", "").replace("", np.nan)
    )

    vehicles_gho = _clean_data(vehicles_gho)

    url = "https://en.wikipedia.org/wiki/List_of_countries_and_territories_by_motor_vehicles_per_capita"
    try:
        vehicles_wiki = pd.read_html(
            url, storage_options=storage_options, encoding="utf8"
        )[0]
        print("File read successfully.")
    except Exception as e:
        logger.warning("Failed to read the file.",e)
        vehicles_wiki = pd.DataFrame(columns=["Country", "country", "number cars"])

    vehicles_wiki.rename(
        columns={"Location": "Country", "Vehicles": "number cars"}, inplace=True
    )

    vehicles_wiki = add_iso2_country_code(vehicles_wiki)

    vehicles_wiki = _clean_data(vehicles_wiki)

    
    # Add missing countries, which are available in the wikipedia source.
    missing_countries = set(vehicles_wiki["country"]) - set(vehicles_gho["country"])
    print(f"Adding the missing countries {missing_countries} from Wikipedia source.")

    vehicles_wiki_to_add = vehicles_wiki[vehicles_wiki["country"].isin(missing_countries)]

    nbr_vehicles = pd.concat([vehicles_gho, vehicles_wiki_to_add], 
        ignore_index=True
    )

    return nbr_vehicles


def download_CO2_emissions():
    """
    Downloads the CO2 emissions from transport in % of total fuel combustion.
    The data is used to estimate the average fuel consumption of land transport.
    It is until the year 2014. # TODO: Maybe search for more recent years or another proxy to 
    estimating the average fuel efficiency (MWh/100km).
    
    The live API of the World Bank has stopped providing the dataset since October 2024. 
    So this link is used: https://web.archive.org/web/20240527231108/https://data.worldbank.org/indicator/EN.CO2.TRAN.ZS?view=map
    """
    url = "https://web.archive.org/web/20240521093243if_/https://api.worldbank.org/v2/en/indicator/EN.CO2.TRAN.ZS?downloadformat=excel"

   # Read the 'Data' sheet directly from the Excel file at the provided URL
    try:
        CO2_emissions = pd.read_excel(url, sheet_name="Data", skiprows=[0, 1, 2])
        print("File read successfully.")
    except Exception as e:
        logger.warning("Failed to read the file. Falling back on hard-coded data:",e)
        return pd.DataFrame()

    CO2_emissions = CO2_emissions[
        ["Country Name", "Country Code", "Indicator Name", "2014"]
    ]

    # Estimate efficiency based on CO2 emissions from transport (% of total fuel combustion)
    CO2_emissions["average fuel efficiency"] = (100 - CO2_emissions["2014"]) / 100

    CO2_emissions = CO2_emissions.rename(columns={"Country Name": "Country"})

    CO2_emissions = add_iso2_country_code(CO2_emissions) 

    return CO2_emissions


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_transport_data_input")

    # configure_logging(snakemake)

    # run = snakemake.config.get("run", {})
    # RDIR = run["name"] + "/" if run.get("name") else ""
    # store_path_data = Path.joinpath(Path().cwd(), "data")
    # country_list = country_list_to_geofk(snakemake.config["countries"])'

    nbr_vehicles = download_number_of_vehicles().copy()

    CO2_emissions = download_CO2_emissions().copy()

    if nbr_vehicles.empty or CO2_emissions.empty:
        # In case one of the urls is not working, we can use the hard-coded data
        src = BASE_DIR + "/data/temp_hard_coded/transport_data.csv"
        dest = snakemake.output.transport_data_input
        shutil.copy(src, dest)
    else:
        # Join the DataFrames by the 'country' column to prepare the tabular transport_data, 
        # which will be saved as transport_data.csv in the resource folder and used
        # to prepare further (nodal) transport data in prepare_transport_data 
        # and to scale the e-mob parameters in prepare_sector_network.
        transport = pd.merge(nbr_vehicles, CO2_emissions, on="country")
        transport = transport[["country", "number cars", "average fuel efficiency"]]

        missing = transport.index[transport["average fuel efficiency"].isna()]
        if not missing.empty:
            print(
                "Missing data on fuel efficiency from:\n",
                f"{list(transport.loc[missing].country)}.",
                "\nFilling gaps with averaged data."
            )

            fill_value = transport["average fuel efficiency"].mean()
            transport.loc[missing, "average fuel efficiency"] = fill_value

        transport.loc[:, "average fuel efficiency"] = transport[
            "average fuel efficiency"
        ].round(3)

        transport.to_csv(
            snakemake.output.transport_data_input,
            sep=",",
            encoding="utf-8",
            header="true",
            index=False,
        )
