# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Retrieve the gas network datasets (GGIT and IGGIELGN) used by
prepare_gas_network.py.
"""

import os
import zipfile
from pathlib import Path

import pandas as pd
from _helpers import BASE_DIR, content_retrieve, create_logger, progress_retrieve

logger = create_logger(__name__)


def download_IGGIELGN_gas_network(fn) -> None:
    """
    Downloads a global dataset for gas networks as .xlsx.

    Data on the European gas transmission network ("gas_network_iggielgn"
    dataset).
    """

    # Save locations
    zip_fn = Path(os.path.join(BASE_DIR, "IGGIELGN.zip"))
    to_fn = Path(os.path.join(BASE_DIR, "data/gas_network/scigrid-gas"))

    logger.info(f"Downloading databundle from '{fn}'.")
    progress_retrieve(fn, zip_fn)

    logger.info(f"Extracting databundle.")
    zipfile.ZipFile(zip_fn).extractall(to_fn)

    zip_fn.unlink()

    logger.info(f"Gas infrastructure data available in '{to_fn}'.")


def download_GGIT_gas_network(fn) -> pd.DataFrame:
    """
    Downloads a global dataset for gas networks as .xlsx.

    Data on gas pipelines worldwide ("pipelines_gem" dataset).
    """
    GGIT_gas_pipeline = pd.read_excel(
        content_retrieve(fn),
        index_col=0,
        sheet_name="Gas Pipelines 2022-12-16",
        header=0,
    )

    return GGIT_gas_pipeline


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_gas_network")

    GGIT_gas_pipeline = download_GGIT_gas_network(snakemake.params.url_ggit)
    GGIT_gas_pipeline.to_csv(snakemake.output.ggit_raw)

    download_IGGIELGN_gas_network(snakemake.params.url_iggielgn)
