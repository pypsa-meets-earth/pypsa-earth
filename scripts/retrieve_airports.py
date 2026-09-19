# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

from _helpers import read_csv_nafix


def download_airports(fn_airports, fn_runways):
    """
    Downloads the world airports as .csv File in addition to runnways
    information.

    Data on airports worldwide ("airports" and "air_runways" datasets).
    """
    storage_options = {"User-Agent": "Mozilla/5.0"}
    airports_csv = read_csv_nafix(
        fn_airports, index_col=0, storage_options=storage_options, encoding="utf8"
    )

    storage_options = {"User-Agent": "Mozilla/5.0"}
    runways_csv = read_csv_nafix(
        fn_runways, index_col=0, storage_options=storage_options, encoding="utf8"
    )

    return (airports_csv, runways_csv)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_airports")

    airports_csv, runways_csv = download_airports(
        snakemake.params.url_airports,
        snakemake.params.url_runways,
    )
    airports_csv.to_csv(snakemake.output.airports_raw)
    runways_csv.to_csv(snakemake.output.runways_raw)
