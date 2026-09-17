# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Prepares airport location and size data used for the aviation sector.

Relevant Settings
-----------------

```yaml
custom_data:
    airports:

sector:
    airport_sizing_factor:
```

Outputs
-------

- ``resources/airports.csv``: medium and large scheduled airports per country. The
  ``fraction`` column gives each airport's share of its country's aviation demand,
  derived from the size weights set by ``airport_sizing_factor``.

Description
-----------

If ``custom_data.airports`` is enabled, the user-provided
``data/custom/airports.csv`` is copied to the output path. Otherwise, the
global OurAirports dataset is downloaded, filtered to commercial airports
with a scheduled service, and combined with runway information. Each
airport is assigned a size weight of 1 if medium and
``sector.airport_sizing_factor`` if large. The ``fraction`` column is that
weight divided by the sum of the weights of all airports in the same
country, so the fractions within a country sum to one. Downstream,
``prepare_sector_network`` multiplies this fraction by the national aviation
demand to obtain the demand served by each airport.

References
----------

- OurAirports, maintained by David Megginson: global airport and runway database,
  released to the public domain (https://ourairports.com/data/). Downloaded from
  the mirror at https://davidmegginson.github.io/ourairports-data/.
"""

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from _helpers import BASE_DIR, read_csv_nafix

# from _helpers import configure_logging


# logger = logging.getLogger(__name__)


def download_airports() -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Download the global airport and runway tables from OurAirports.

    The tables are fetched live from the OurAirports mirror, so the number of
    airports and runways they contain changes as the upstream dataset is
    updated. See the module-level References section for the data source.

    Returns
    -------
    airports_csv : pd.DataFrame
        World airports, one row per airport.
    runways_csv : pd.DataFrame
        World runways, one row per runway.
    """
    fn = "https://davidmegginson.github.io/ourairports-data/airports.csv"
    storage_options = {"User-Agent": "Mozilla/5.0"}
    airports_csv = read_csv_nafix(
        fn, index_col=0, storage_options=storage_options, encoding="utf8"
    )

    fn = "https://davidmegginson.github.io/ourairports-data/runways.csv"
    storage_options = {"User-Agent": "Mozilla/5.0"}
    runways_csv = read_csv_nafix(
        fn, index_col=0, storage_options=storage_options, encoding="utf8"
    )

    return (airports_csv, runways_csv)


def preprocess_airports(df: pd.DataFrame) -> pd.DataFrame:
    """
    Preprocess the airports data.

    Parameters
    ----------
    df : pd.DataFrame
        Merged airports and runways data, as returned by :func:`download_airports`.

    Returns
    -------
    pd.DataFrame
        Medium and large scheduled, commercial airports. The ``fraction``
        column holds each airport's size weight divided by the total weight
        of all airports in the same country, and ``iso_country`` is renamed
        to ``country``.
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
