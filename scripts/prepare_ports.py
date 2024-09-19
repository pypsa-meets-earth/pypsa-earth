# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
import logging
import os
from pathlib import Path

import country_converter as coco
import numpy as np
import pandas as pd

# from _helpers import configure_logging


# logger = logging.getLogger(__name__)


def download_ports():
    """
    Downloads the world ports index csv File and NOT as shape or other because
    it is updated on a monthly basis.

    The following csv file was downloaded from the webpage
    https://msi.nga.mil/Publications/WPI
    as a csv file that is updated monthly as mentioned on the webpage. The dataset contains 3711 ports.
    """
    fn = "https://msi.nga.mil/api/publications/download?type=view&key=16920959/SFH00000/UpdatedPub150.csv"
    wpi_csv = pd.read_csv(fn, index_col=0)

    return wpi_csv


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_ports")

    # configure_logging(snakemake)

    # run = snakemake.config.get("run", {})
    # RDIR = run["name"] + "/" if run.get("name") else ""
    # store_path_data = Path.joinpath(Path().cwd(), "data")
    # country_list = country_list_to_geofk(snakemake.config["countries"])'

    df = download_ports().copy()

    # Add ISO2 country code for each country
    df = df.rename(
        columns={
            "Country Code": "country_full_name",
            "Latitude": "y",
            "Longitude": "x",
            "Main Port Name": "name",
        }
    )
    df["country"] = df.country_full_name.apply(
        lambda x: coco.convert(names=x, to="ISO2", not_found=None)
    )

    # Drop small islands that have no ISO2:
    df = df[df.country_full_name != "Wake Island"]
    df = df[df.country_full_name != "Johnson Atoll"]
    df = df[df.country_full_name != "Midway Islands"]

    # Select the columns that we need to keep
    df = df.reset_index()
    df = df[
        [
            "World Port Index Number",
            "Region Name",
            "name",
            "Alternate Port Name",
            "country",
            "World Water Body",
            "Liquified Natural Gas Terminal Depth (m)",
            "Harbor Size",
            "Harbor Type",
            "Harbor Use",
            "country_full_name",
            "y",
            "x",
        ]
    ]

    # Drop ports that are very small and that have unknown size (Unknown size ports are in total 19 and not suitable for H2 - checked visually)
    ports = df.loc[df["Harbor Size"].isin(["Small", "Large", "Medium"])]

    ports.insert(8, "Harbor_size_nr", 1)
    ports.loc[ports["Harbor Size"].isin(["Small"]), "Harbor_size_nr"] = 1
    ports.loc[ports["Harbor Size"].isin(["Medium"]), "Harbor_size_nr"] = 2
    ports.loc[ports["Harbor Size"].isin(["Large"]), "Harbor_size_nr"] = 3

    df1 = ports.copy()
    df1 = df1.groupby(["country_full_name"]).sum("Harbor_size_nr")
    df1 = df1[["Harbor_size_nr"]]
    df1 = df1.rename(columns={"Harbor_size_nr": "Total_Harbor_size_nr"})

    ports = ports.set_index("country_full_name").join(df1, how="left")

    ports["fraction"] = ports["Harbor_size_nr"] / ports["Total_Harbor_size_nr"]

    ports.to_csv(snakemake.output[0], sep=",", encoding="utf-8", header="true")
