# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

from _helpers import read_csv_nafix


def download_ports(fn):
    """
    Downloads the world ports index csv File and NOT as shape or other because
    it is updated on a monthly basis.

    Data on sea ports worldwide ("sea_ports_nga" dataset).
    """
    wpi_csv = read_csv_nafix(fn, index_col=0)

    return wpi_csv


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_ports")

    wpi_csv = download_ports(snakemake.params.url_ports)
    wpi_csv.to_csv(snakemake.output.ports_raw)
