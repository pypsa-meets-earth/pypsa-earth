# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Build historical annual ammonia production per country in ktonNH3/a.

Description
-------

This functions takes data from the `Minerals Yearbook <https://www.usgs.gov/centers/national-minerals-information-center/nitrogen-statistics-and-information>`_
 (July 2024) published by the US Geological Survey (USGS) and the National Minerals Information Center and extracts the annual ammonia production per country in ktonN/a. The data is converted to ktonNH3/a.
"""

import country_converter as cc
import pandas as pd
from _helpers import configure_logging, content_retrieve, create_logger

logger = create_logger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_ammonia_production")

    configure_logging(snakemake)

    # Download ammonia production data - try primary source, fallback to mirror if needed
    logger.info("Downloading ammonia production data from USGS...")
    try:
        url = "https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/media/files/myb1-2023-nitro-ERT.xlsx"
        content = content_retrieve(url)
    except Exception:
        logger.warning("Primary source failed, trying fallback mirror...")
        url = "https://data.pypsa.org/workflows/eur/nitrogen_statistics/2023/myb1-2023-nitro-ERT.xlsx"
        content = content_retrieve(url)

    ammonia = pd.read_excel(
        content,
        sheet_name="T12",
        skiprows=5,
        header=0,
        index_col=0,
        skipfooter=7,
        na_values=["--"],
    )

    ammonia.index = cc.convert(ammonia.index, to="iso2")

    years = [str(i) for i in range(2019, 2024)]

    ammonia = ammonia.rename(columns=lambda x: str(x))[years]

    # convert from ktonN to ktonNH3
    ammonia *= 17 / 14

    ammonia.index.name = "ktonNH3/a"

    ammonia.to_csv(snakemake.output.ammonia_production)
