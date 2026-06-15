# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Retrieve EDGAR CO2 emission data from the JRC EDGAR database.

Downloads the EDGAR v8.0 CO2 fossil fuel emission dataset (excl. short-cycle
organic carbon) and saves it as a local XLSX file for subsequent processing
by build_co2_emissions.
"""

from _helpers import configure_logging, content_retrieve, create_logger

logger = create_logger(__name__)

EDGAR_CO2_SOURCES = {
    "primary": "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/v80_GHG/CO2_excl_short-cycle_org_C/EDGAR_v80_GHG_CO2_excl_short-cycle_org_C_1970_2022.xlsx",
    "archive": "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/v80_FT2022_GHG/CO2_excl_short-cycle_org_C/EDGAR_v80_FT2022_GHG_CO2_excl_short-cycle_org_C_1970_2022.xlsx",
}


def download_edgar_co2_data(sources: dict) -> bytes:
    """
    Download EDGAR CO2 fossil emission data, with fallback URL on failure.

    Parameters
    ----------
    sources : dict
        Dictionary with 'primary' and 'archive' URL keys.

    Returns
    -------
    io.BytesIO
        The content of the downloaded EDGAR XLSX file.
    """
    logger.info("Downloading EDGAR CO2 emission data (primary source)...")
    try:
        content = content_retrieve(sources["primary"])
        logger.info("Downloaded EDGAR CO2 data from primary source.")
        return content
    except Exception as e:
        logger.warning(f"Primary source failed ({e}), trying fallback...")
        content = content_retrieve(sources["archive"])
        logger.info("Downloaded EDGAR CO2 data from fallback source.")
        return content


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_emissions")

    configure_logging(snakemake)

    content = download_edgar_co2_data(EDGAR_CO2_SOURCES)

    with open(snakemake.output.edgar, "wb") as f:
        f.write(content.read())

    logger.info(f"EDGAR CO2 data saved to '{snakemake.output.edgar}'.")
