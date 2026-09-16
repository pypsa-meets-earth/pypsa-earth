# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Retrieve the number-of-vehicles and CO2-emissions-from-transport datasets
used by prepare_transport_data_input.py.
"""

import logging

import country_converter as coco
import numpy as np
import pandas as pd
from _helpers import read_csv_nafix

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
    df = df.loc[df.country.apply(lambda x: isinstance(x, str)), :]

    return df


def download_number_of_vehicles(fn_who, fn_wiki):
    """
    Downloads and returns the number of registered vehicles as tabular data
    from the Global Health Observatory (GHO) repository data and from Wikipedia.

    Data on the number of registered motor vehicles worldwide
    ("n_vehicles_who" and "vehicles_per_capita_wiki" datasets).
    """

    def _clean_data(df):
        df = df.dropna(subset=["number cars"])
        df.loc[:, "number cars"] = df.loc[:, "number cars"].astype(int)
        return df  # [["Country", "number cars"]]

    storage_options = {"User-Agent": "Mozilla/5.0"}
    try:
        vehicles_gho = read_csv_nafix(
            fn_who, storage_options=storage_options, encoding="utf8"
        )
        logger.info("File read successfully.")
    except Exception as e:
        logger.warning(
            f"Failed to read the file. Falling back on hard-coded data. \nError: {e}"
        )
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

    try:
        vehicles_wiki = pd.read_html(
            fn_wiki, storage_options=storage_options, encoding="utf8"
        )[0]
        logger.info("File read successfully.")
    except Exception as e:
        logger.warning("Failed to read the file.", e)
        vehicles_wiki = pd.DataFrame(columns=["Country", "country", "number cars"])

    vehicles_wiki.rename(
        columns={"Region": "Country", "Road motor vehicles": "number cars"},
        inplace=True,
    )

    vehicles_wiki = add_iso2_country_code(vehicles_wiki)

    vehicles_wiki = _clean_data(vehicles_wiki)

    # Add missing countries, which are available in the wikipedia source.
    missing_countries = set(vehicles_wiki["country"]) - set(vehicles_gho["country"])
    logger.info(
        f"Adding the missing countries {missing_countries} from Wikipedia source."
    )

    vehicles_wiki_to_add = vehicles_wiki[
        vehicles_wiki["country"].isin(missing_countries)
    ]

    nbr_vehicles = pd.concat([vehicles_gho, vehicles_wiki_to_add], ignore_index=True)

    return nbr_vehicles


def download_CO2_emissions(fn):
    """
    Downloads the CO2 emissions from transport in % of total fuel combustion.
    The data is used to estimate the average fuel consumption of land transport.
    It is until the year 2014. # TODO: Maybe search for more recent years or another proxy to
    estimating the average fuel efficiency (MWh/100km).

    Data on CO2 emissions from transport worldwide
    ("transport_emission_worldbank" dataset).
    """
    # Read the 'Data' sheet directly from the Excel file at the provided URL
    try:
        CO2_emissions = pd.read_excel(fn, sheet_name="Data", skiprows=[0, 1, 2])
        logger.info("File read successfully.")
    except Exception as e:
        logger.warning("Failed to read the file. Falling back on hard-coded data:", e)
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

        snakemake = mock_snakemake("retrieve_transport_data_input")

    nbr_vehicles = download_number_of_vehicles(
        snakemake.params.url_n_vehicles_who,
        snakemake.params.url_vehicles_per_capita_wiki,
    )
    nbr_vehicles.to_csv(snakemake.output.n_vehicles_raw, index=False)

    CO2_emissions = download_CO2_emissions(
        snakemake.params.url_transport_emission_worldbank
    )
    CO2_emissions.to_csv(snakemake.output.transport_emissions_raw, index=False)
