# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from _helpers import BASE_DIR

# from _helpers import configure_logging


# logger = logging.getLogger(__name__)


def download_airports():
    """
    Downloads the world airports as .csv File in addition to runnways
    information.

    The following csv file was downloaded from the webpage
    https://ourairports.com/data/
    as a .csv file. The dataset contains 74844 airports.
    """
    fn = "https://davidmegginson.github.io/ourairports-data/airports.csv"
    storage_options = {"User-Agent": "Mozilla/5.0"}
    airports_csv = pd.read_csv(
        fn, index_col=0, storage_options=storage_options, encoding="utf8"
    )

    fn = "https://davidmegginson.github.io/ourairports-data/runways.csv"
    storage_options = {"User-Agent": "Mozilla/5.0"}
    runways_csv = pd.read_csv(
        fn, index_col=0, storage_options=storage_options, encoding="utf8"
    )

    return (airports_csv, runways_csv)


def preprocess_airports(df):
    """
    Preprocess the airports data
    """

    # Keep only airports that are of type medium and large
    df = df.loc[df["type"].isin(["large_airport", "medium_airport"])]

    # Filtering out the military airbases and keeping only commercial airports
    df = df[~df.iata_code.isnull()]

    # Keep only airports that have schedules
    df = df.loc[df["scheduled_service"].isin(["yes"])]

    df.insert(2, "airport_size_nr", 1)
    df.loc[df["type"].isin(["medium_airport"]), "airport_size_nr"] = 1
    df.loc[df["type"].isin(["large_airport"]), "airport_size_nr"] = (
        snakemake.params.airport_sizing_factor
    )

    # Calculate the number of total airports size
    df1 = df.copy()
    df1 = df1.groupby(["iso_country"]).sum("airport_size_nr")
    df1 = df1[["airport_size_nr"]]
    df1 = df1.rename(columns={"airport_size_nr": "Total_airport_size_nr"}).reset_index()

    # Merge dataframes to get additional info on runnway for most ports
    airports = pd.merge(
        df, df1, how="left", left_on="iso_country", right_on="iso_country"
    )

    # Calculate fraction based on size
    airports["fraction"] = (
        airports["airport_size_nr"] / airports["Total_airport_size_nr"]
    )

    # Rename columns
    airports = airports.rename(columns={"iso_country": "country"})

    return airports


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_airports")
    # configure_logging(snakemake)

    # run = snakemake.config.get("run", {})
    # RDIR = run["name"] + "/" if run.get("name") else ""
    # store_path_data = Path.joinpath(Path().cwd(), "data")
    # country_list = country_list_to_geofk(snakemake.config["countries"])'

    if snakemake.params.airport_custom_data:
        custom_airports = Path(BASE_DIR).joinpath("data", "custom", "airports.csv")
        shutil.copy(custom_airports, snakemake.output[0])
    else:
        # Prepare downloaded data
        download_data = download_airports()

        airports_csv = download_data[0].copy()
        airports_csv = airports_csv[
            [
                "ident",
                "type",
                "name",
                "latitude_deg",
                "longitude_deg",
                "elevation_ft",
                "continent",
                "iso_country",
                "iso_region",
                "municipality",
                "scheduled_service",
                "iata_code",
            ]
        ]
        airports_csv.loc[airports_csv["iso_country"].isnull(), "iso_country"] = "NA"
        airports_csv = airports_csv.rename(columns={"latitude_deg": "y"})
        airports_csv = airports_csv.rename(columns={"longitude_deg": "x"})

        runways_csv = download_data[1].copy()
        runways_csv = runways_csv[
            ["airport_ident", "length_ft", "width_ft", "surface", "lighted", "closed"]
        ]
        runways_csv = runways_csv.drop_duplicates(subset=["airport_ident"])

        airports_original = pd.merge(
            airports_csv,
            runways_csv,
            how="left",
            left_on="ident",
            right_on="airport_ident",
        )
        airports_original = airports_original.drop("airport_ident", axis=1)

        df = airports_original.copy()

        airports = preprocess_airports(df)

        # Save
        airports.to_csv(snakemake.output[0], sep=",", encoding="utf-8", header="true")
