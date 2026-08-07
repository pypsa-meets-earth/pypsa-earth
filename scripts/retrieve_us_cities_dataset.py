# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Retrieve US cities dataset.
"""

import zipfile

import pandas as pd
from _helpers import configure_logging, content_retrieve, create_logger, to_csv_nafix

logger = create_logger(__name__)

US_CITIES_SOURCE = "https://simplemaps.com/static/data/us-cities/1.93/basic/simplemaps_uscities_basicv1.93.zip"


def load_us_cities_data(us_cities_source: str) -> pd.DataFrame:
    """
    Load US cities data from a ZIP file containing CSV.

    Parameters
    ----------
    us_cities_source : str
        URL to the US cities ZIP file.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing US cities with columns: city, state_id, lat, lng, population, etc.
    """
    # Download the ZIP file
    zip_content = content_retrieve(us_cities_source)

    # Extract and read the CSV from the ZIP
    with zipfile.ZipFile(zip_content) as z:
        # Get the CSV filename (usually 'uscities.csv' or similar)
        csv_filename = [f for f in z.namelist() if f.endswith(".csv")][0]

        with z.open(csv_filename) as csv_file:
            us_cities = pd.read_csv(csv_file, encoding="utf-8")

    logger.info(f"Loaded {len(us_cities)} US cities")

    # Lowercase city names for matching
    us_cities["city"] = us_cities["city"].str.lower()

    return us_cities


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_us_cities_dataset")

    configure_logging(snakemake)

    # Load US cities data
    us_cities = load_us_cities_data(US_CITIES_SOURCE)

    # Save the results
    to_csv_nafix(us_cities, snakemake.output.us_cities, index=False)
