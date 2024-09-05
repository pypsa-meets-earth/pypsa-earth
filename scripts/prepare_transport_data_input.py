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

# from _helpers import configure_logging


# logger = logging.getLogger(__name__)


def download_number_of_vehicles():
    """
    Downloads the Number of registered vehicles as .csv File.

    The following csv file was downloaded from the webpage
    https://apps.who.int/gho/data/node.main.A995
    as a .csv file.
    """
    fn = "https://apps.who.int/gho/athena/data/GHO/RS_194?filter=COUNTRY:*&ead=&x-sideaxis=COUNTRY;YEAR;DATASOURCE&x-topaxis=GHO&profile=crosstable&format=csv"
    storage_options = {"User-Agent": "Mozilla/5.0"}

    # Read the 'Data' sheet directly from the csv file at the provided URL
    try:
        Nbr_vehicles_csv = pd.read_csv(
            fn, storage_options=storage_options, encoding="utf8"
        )
        print("File read successfully.")
    except Exception as e:
        print("Failed to read the file:", e)
        return pd.DataFrame()

    Nbr_vehicles_csv = Nbr_vehicles_csv.rename(
        columns={
            "Countries, territories and areas": "Country",
            "Number of registered vehicles": "number cars",
        }
    )

    # Add ISO2 country code for each country
    cc = coco.CountryConverter()
    Country = pd.Series(Nbr_vehicles_csv["Country"])
    Nbr_vehicles_csv["country"] = cc.pandas_convert(
        series=Country, to="ISO2", not_found="not found"
    )

    # # Remove spaces, Replace empty values with NaN
    Nbr_vehicles_csv["number cars"] = (
        Nbr_vehicles_csv["number cars"].str.replace(" ", "").replace("", np.nan)
    )

    # Drop rows with NaN values in 'number cars'
    Nbr_vehicles_csv = Nbr_vehicles_csv.dropna(subset=["number cars"])

    # convert the 'number cars' to integer
    Nbr_vehicles_csv["number cars"] = Nbr_vehicles_csv["number cars"].astype(int)

    return Nbr_vehicles_csv


def download_CO2_emissions():
    """
    Downloads the CO2_emissions from vehicles as .csv File.

    The dataset is downloaded from the following link: https://data.worldbank.org/indicator/EN.CO2.TRAN.ZS?view=map
    It is until the year 2014. # TODO: Maybe search for more recent years.
    """
    url = (
        "https://api.worldbank.org/v2/en/indicator/EN.CO2.TRAN.ZS?downloadformat=excel"
    )

    # Read the 'Data' sheet directly from the Excel file at the provided URL
    try:
        CO2_emissions = pd.read_excel(url, sheet_name="Data", skiprows=[0, 1, 2])
        print("File read successfully.")
    except Exception as e:
        print("Failed to read the file:", e)
        return pd.DataFrame()

    CO2_emissions = CO2_emissions[
        ["Country Name", "Country Code", "Indicator Name", "2014"]
    ]

    # Calculate efficiency based on CO2 emissions from transport (% of total fuel combustion)
    CO2_emissions["average fuel efficiency"] = (100 - CO2_emissions["2014"]) / 100

    # Add ISO2 country code for each country
    CO2_emissions = CO2_emissions.rename(columns={"Country Name": "Country"})
    cc = coco.CountryConverter()
    Country = pd.Series(CO2_emissions["Country"])
    CO2_emissions["country"] = cc.pandas_convert(
        series=Country, to="ISO2", not_found="not found"
    )

    # Drop region names that have no ISO2:
    CO2_emissions = CO2_emissions[CO2_emissions.country != "not found"]

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

    # Downloaded and prepare vehicles_csv:
    vehicles_csv = download_number_of_vehicles().copy()

    # Downloaded and prepare CO2_emissions_csv:
    CO2_emissions_csv = download_CO2_emissions().copy()

    if vehicles_csv.empty or CO2_emissions_csv.empty:
        # In case one of the urls is not working, we can use the hard-coded data
        src = os.getcwd() + "/data/temp_hard_coded/transport_data.csv"
        dest = snakemake.output.transport_data_input
        shutil.copy(src, dest)
    else:
        # Join the DataFrames by the 'country' column
        merged_df = pd.merge(vehicles_csv, CO2_emissions_csv, on="country")
        merged_df = merged_df[["country", "number cars", "average fuel efficiency"]]

        # Drop rows with NaN values in 'average fuel efficiency'
        merged_df = merged_df.dropna(subset=["average fuel efficiency"])

        # Convert the 'average fuel efficiency' to float
        merged_df["average fuel efficiency"] = merged_df[
            "average fuel efficiency"
        ].astype(float)

        # Round the 'average fuel efficiency' to three decimal places
        merged_df.loc[:, "average fuel efficiency"] = merged_df[
            "average fuel efficiency"
        ].round(3)

        # Save the merged DataFrame to a CSV file
        merged_df.to_csv(
            snakemake.output.transport_data_input,
            sep=",",
            encoding="utf-8",
            header="true",
            index=False,
        )
