# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Prepares the urban population share per country, used to split population-
dependent demand between urban and rural nodes.

Outputs
-------

- ``resources/urban_percent.csv``: urban population as a percentage of total population, per country (ISO2 code) and year (2020 onwards).

Description
-----------

Downloads the UNCTAD "Total and urban population, annual" dataset, keeps years
from 2020 onwards, maps country/economy names to ISO2 codes, and averages the
urban percentage where multiple entries map to the same country and year.

References
----------

- Total and urban population, annual (``US.PopTotal``), United Nations Trade and
  Development Data Hub (UNCTADstat), https://unctadstat.unctad.org/datacentre/.
  Licensed under Creative Commons Attribution 3.0 IGO, with UNCTAD asking that
  "United Nations Trade and Development Data Hub" be cited as the source
  (https://unctadstat.unctad.org/EN/About.html#TERMS_OF_USE).
"""

import os

import country_converter as coco
import pandas as pd
import py7zr
import requests
from _helpers import read_csv_nafix

# from _helpers import configure_logging


# logger = logging.getLogger(__name__)


def download_urban_percent() -> pd.DataFrame:
    """
    Downloads the UNCTAD "Total and urban population, annual" .7z File and
    extracts it as csv File.

    The dataset is UNCTADstat series ``US.PopTotal``, browsable at
    https://unctadstat.unctad.org/datacentre/ and retrieved here through the
    bulk download endpoint, which serves it as ``US_PopTotal.csv.7z``. It
    covers urban percent for most countries from 1950 with projections to
    2050; this script keeps 2020 onwards. See the module-level References
    section for the citation requested by UNCTAD.

    Returns
    -------
    pd.DataFrame
        Total and urban population per country and year, as published by
        UNCTAD.
    """
    url = "https://unctadstat-api.unctad.org/bulkdownload/US.PopTotal/US_PopTotal"

    # Make a GET request to the URL
    response = requests.get(url)

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        # Extract the filename from the Content-Disposition header
        content_disposition = response.headers.get("Content-Disposition")
        if content_disposition:
            filename = content_disposition.split("filename=")[1].strip('"')
        else:
            filename = "downloaded_file.csv.7z"  # Provide a default filename if Content-Disposition header is not present

        # Write the content of the response to a file
        with open(filename, "wb") as f:
            f.write(response.content)

        print(f"Urban percent downloaded successfully as {filename}")

        # Extract the downloaded .7z file
        with py7zr.SevenZipFile(filename, "r") as archive:
            archive.extractall()

        print(f"Urban percent extracted successfully")

        # Read the extracted CSV file
        csv_filename = os.path.splitext(filename)[
            0
        ]  # Remove the .7z extension to get the CSV filename
        urban_percent_orig = read_csv_nafix(csv_filename)

        print("Urban percent CSV file read successfully:")

        # Remove the downloaded .7z and .csv files
        os.remove(filename)
        os.remove(csv_filename)

    else:
        print(f"Failed to download file: Status code {response.status_code}")

    return urban_percent_orig


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_urban_percent")

    df = download_urban_percent().copy()

    # Select the columns that we need to keep
    df = df[
        [
            "Year",
            "Economy Label",
            "Urban population as percentage of total population",
        ]
    ]

    # Keep only years above 2020
    df = df.loc[(df["Year"] >= 2020)]

    # Add ISO2 country code for each country
    cc = coco.CountryConverter()
    Economy_Label = pd.Series(df["Economy Label"])
    df["country"] = cc.pandas_convert(
        series=Economy_Label, to="ISO2", not_found="not found"
    )

    # Drop isos that were not found:
    df = df.loc[df["country"] != "not found"]

    # Drop region names where country column contains list of countries
    df = df.loc[df.country.apply(lambda x: isinstance(x, str)), :]

    # Reduce the data to one value for the urban percent per country and year
    df = df.groupby(["country", "Year"], as_index=False).mean(numeric_only=True)

    df = df.set_index("country")

    # Save
    df.to_csv(snakemake.output[0], sep=",", encoding="utf-8", header="true")
