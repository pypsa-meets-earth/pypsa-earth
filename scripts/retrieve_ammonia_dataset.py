# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Retrieve ammonia production dataset from USGS.
"""

from _helpers import configure_logging, content_retrieve, create_logger

logger = create_logger(__name__)

USGS_AMMONIA_SOURCES = {
    "primary": "https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/media/files/myb1-2023-nitro-ERT.xlsx",
    "archive": "https://data.pypsa.org/workflows/eur/nitrogen_statistics/2023/myb1-2023-nitro-ERT.xlsx",
}


def download_ammonia_production_data(ammonia_sources: dict) -> bytes:
    """
    Download ammonia production data from the USGS website, with a fallback to an archived version if the primary source is unavailable.
    Parameters
    ----------
    ammonia_sources : dict
        Dictionary containing the URLs for the primary and archive sources.
    Returns
    -------
    bytes
        The content of the downloaded ammonia production data.
    """
    # Download ammonia production data - try primary source, fallback to mirror if needed
    logger.info("Downloading ammonia production data from USGS...")
    try:
        url = ammonia_sources["primary"]
        content = content_retrieve(url)
    except Exception:
        logger.warning("Primary source failed, trying fallback mirror...")
        url = ammonia_sources["archive"]
        content = content_retrieve(url)

    return content


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_ammonia_production")

    configure_logging(snakemake)

    # Download ammonia production dataset
    content = download_ammonia_production_data(USGS_AMMONIA_SOURCES)

    # Save raw Excel file
    with open(snakemake.output.usgs_ammonia_dataset, "wb") as f:
        f.write(content.read())
