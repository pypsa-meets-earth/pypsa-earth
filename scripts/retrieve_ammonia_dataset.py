# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Retrieve ammonia production dataset from USGS.
"""

from _helpers import configure_logging, content_retrieve, create_logger

logger = create_logger(__name__)


def download_ammonia_production_data(url_primary: str, url_archive: str) -> bytes:
    """
    Download ammonia production data from the USGS website, with a fallback to an archived version if the primary source is unavailable.
    Parameters
    ----------
    url_primary : str
        URL of the primary USGS source.
    url_archive : str
        URL of the archived fallback mirror.
    Returns
    -------
    bytes
        The content of the downloaded ammonia production data.
    """
    # Download ammonia production data - try primary source, fallback to mirror if needed
    logger.info("Downloading ammonia production data from USGS...")
    try:
        content = content_retrieve(url_primary)
    except Exception:
        logger.warning("Primary source failed, trying fallback mirror...")
        content = content_retrieve(url_archive)

    return content


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_ammonia_production")

    configure_logging(snakemake)

    # Download ammonia production dataset
    content = download_ammonia_production_data(
        snakemake.params.url_primary,
        snakemake.params.url_archive,
    )

    # Save raw Excel file
    with open(snakemake.output.usgs_ammonia_dataset, "wb") as f:
        f.write(content.read())
