# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Prepares port location and size data used for the shipping and hydrogen
export sectors.

Relevant Settings
-----------------

```yaml
custom_data:
    export_ports:
```

Outputs
-------

- ``resources/ports.csv``: world ports per country. The ``fraction`` column gives
  each port's share of its country's total harbor size weighting.
- ``resources/export_ports.csv``: the ports of each country's largest available
  harbor size class, used as candidate hydrogen/derivative export locations; or the
  user-provided ``data/custom/export_ports.csv`` if ``custom_data.export_ports`` is
  enabled.

Description
-----------

Downloads the World Port Index, drops entries that could not be matched to
an ISO2 country code (e.g. small islands), and weights each port by harbor
size as 3 (large), 2 (medium) or 1 (small). The ``fraction`` column is that
weight divided by the sum of the weights of all ports in the same country, so
the fractions within a country sum to one. Candidate export ports are then
selected per country by taking all ports of the largest harbor size present in
that country: all large ports if any exist, otherwise all medium ports,
otherwise all small ones.

References
----------

- World Port Index (Publication 150), U.S. National Geospatial-Intelligence
  Agency, Maritime Safety Office, updated monthly
  (https://msi.nga.mil/Publications/WPI). NGA states that information it
  presents is public information and may be distributed or copied unless
  otherwise specified, with appropriate credit requested
  (https://www.nga.mil/resources/Privacy_Policy.html).
"""

import logging
import os
import shutil
from pathlib import Path

import country_converter as coco
import numpy as np
import pandas as pd
from _helpers import BASE_DIR, read_csv_nafix

# from _helpers import configure_logging


# logger = logging.getLogger(__name__)


def download_ports() -> pd.DataFrame:
    """
    Download the World Port Index as a csv file.

    The csv format is used rather than a shapefile or other format because the
    publication is updated monthly. The table is fetched live, so the number of
    ports it contains changes as the upstream publication is updated. See the
    module-level References section for the data source.

    Returns
    -------
    pd.DataFrame
        World Port Index, one row per port.
    """
    fn = "https://msi.nga.mil/api/publications/download?type=view&key=16920959/SFH00000/UpdatedPub150.csv"
    wpi_csv = read_csv_nafix(fn, index_col=0)

    return wpi_csv


def filter_ports(dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Select, for each country, the ports of the largest harbor size present in
    that country.

    Countries with at least one large port contribute all of their large ports;
    countries without any large port contribute all of their medium ports, and
    countries with neither contribute all of their small ports. A country may
    therefore appear several times in the result.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Ports data with a 'Harbor Size' column ('Large', 'Medium' or 'Small')
        and a 'country' column.

    Returns
    -------
    pd.DataFrame
        All ports of each country's largest available harbor size class, with
        one row per port.
    """
    # Filter large sized ports
    large_ports = dataframe[dataframe["Harbor Size"] == "Large"]
    countries_with_large_ports = large_ports["country"].unique()

    # Filter out countries with large ports
    remaining_ports = dataframe[~dataframe["country"].isin(countries_with_large_ports)]

    # Filter medium sized ports from remaining ports
    medium_ports = remaining_ports[remaining_ports["Harbor Size"] == "Medium"]
    countries_with_medium_ports = medium_ports["country"].unique()

    # Filter out countries with medium ports
    remaining_ports = remaining_ports[
        ~remaining_ports["country"].isin(countries_with_medium_ports)
    ]

    # Filter small sized ports from remaining ports
    small_ports = remaining_ports[remaining_ports["Harbor Size"] == "Small"]

    # Combine all filtered ports
    filtered_ports = pd.concat([large_ports, medium_ports, small_ports])

    return filtered_ports


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_ports")

    config = snakemake.config
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

    if snakemake.params.custom_export:
        custom_export_path = Path(BASE_DIR).joinpath(
            "data", "custom", "export_ports.csv"
        )
        shutil.copy(custom_export_path, snakemake.output[1])

    else:
        filter_ports(ports).to_csv(
            snakemake.output[1], sep=",", encoding="utf-8", header="true"
        )
