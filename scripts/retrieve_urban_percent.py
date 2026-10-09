# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import os

import py7zr
import requests
from _helpers import read_csv_nafix


def download_urban_percent(fn):
    """
    Downloads the United Nations "Total and urban population, annual" .7z File
    and extracts it as csv File.

    The above file was downloaded from the webpage
    https://unctadstat.unctad.org/datacentre/
    as a .7z file. The dataset contains urban percent for most countries from 1950 and predictions until 2050.
    """
    # Make a GET request to the URL
    response = requests.get(fn)

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

        snakemake = mock_snakemake("retrieve_urban_percent")

    df = download_urban_percent(snakemake.params.url_urban_percent)
    df.to_csv(snakemake.output.urban_percent_raw)
