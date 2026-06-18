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
import country_converter as coco
import pandas as pd
from _helpers import configure_logging, create_logger, read_csv_nafix, to_csv_nafix

logger = create_logger(__name__)

# City corrections for known data quality issues in USGS dataset
CITY_CORRECTIONS = {
    "faustina (donaldsonville)": "donaldsonville",
    "greenville": "greeneville",
    "pryor": "pryor creek",
    "geismar": "gonzales",  # closest town to geismar
    "port neal": "sergeant bluff",  # closest town to port neal
}


def clean_city_name(location: str) -> str | None:
    """
    Extract and clean city name from location string.
    Applies manual corrections for known data quality issues.

    Parameters
    ----------
    location : str
        Location string in format "City, State" (e.g., "Donaldsonville, LA").

    Returns
    -------
    str | None
        Cleaned lowercase city name, or None if location is null.
    """
    if pd.isnull(location):
        return None

    # Extract city part (before comma)
    city = location.split(", ")[0].lower().strip()

    # Apply manual corrections
    city = CITY_CORRECTIONS.get(city, city)

    return city


def extract_state(location: str) -> str | None:
    """
    Extract state code from location string.

    Parameters
    ----------
    location : str
        Location string in format "City, State" (e.g., "Donaldsonville, Louisiana").

    Returns
    -------
    str | None
        Two-letter uppercase state code (e.g., "LA"), or None if location is null.
    """
    if pd.isnull(location):
        return None

    parts = location.split(", ")
    if len(parts) > 1:
        return parts[1][:2].upper()
    return None


def extract_country_level_production(content: bytes) -> pd.DataFrame:
    """
    Extract country-level ammonia production totals from the downloaded content.
    Parameters
    ----------
    content : bytes
        The content of the downloaded ammonia production data.
    Returns
    -------
    pd.DataFrame
        A DataFrame containing the country-level ammonia production totals.
    """
    # Extract country-level production totals
    ammonia = pd.read_excel(
        content,
        sheet_name="T12",
        skiprows=5,
        header=0,
        index_col=0,
        skipfooter=7,
        na_values=["--"],
    )

    # Convert country names to ISO2 codes
    cc = coco.CountryConverter()
    ammonia.index = cc.convert(ammonia.index, to="iso2")

    years = [str(i) for i in range(2019, 2024)]
    ammonia = ammonia.rename(columns=lambda x: str(x))[years]

    # Convert from ktonN to ktonNH3
    ammonia *= 17 / 14
    ammonia.index.name = "ktonNH3/a"

    return ammonia


def extract_us_plant_data(content: bytes, us_cities: pd.DataFrame) -> pd.DataFrame:
    """
    Extract US ammonia plant data from the downloaded content.
    Parameters
    ----------
    content : bytes
        The content of the downloaded ammonia production data.
    us_cities : pd.DataFrame
        A DataFrame containing US cities data.
    Returns
    -------
    pd.DataFrame
        A DataFrame containing the US ammonia plant data.
    """
    # Extract US plant data
    us_plants = pd.read_excel(
        content,
        sheet_name="T4",
        skiprows=5,  # Adjust based on actual structure
        header=0,
        na_values=["--"],
        skipfooter=5,
    )

    # Handle "Do." (ditto) notation - replace with actual company name from previous row
    us_plants["Company"] = us_plants["Company"].replace("Do.", pd.NA).ffill()

    # Get city and state from location using helper functions
    us_plants["city"] = us_plants["Location"].apply(clean_city_name)
    us_plants["state"] = us_plants["Location"].apply(extract_state)

    # Merge with US cities data to get lat/lon
    us_plants = us_plants.merge(
        us_cities, how="left", left_on=["city", "state"], right_on=["city", "state_id"]
    )

    # Set plant name as combination of company and city for uniqueness
    us_plants["plant"] = us_plants["Company"] + " - " + us_plants["city"].str.title()

    # Get relevant columns and rename for consistency
    us_plants = us_plants[["plant", "Capacity2", "lat", "lng"]].rename(
        columns={"Capacity2": "capacity", "lat": "y", "lng": "x"}
    )

    us_plants.index.name = "ktonNH3/a"
    us_plants["country"] = "US"

    return us_plants


def extract_plant_data(ammonia_plants_path: str) -> pd.DataFrame:
    """
    Extract ammonia plant data from a provided CSV file.

    Parameters
    ----------
    ammonia_plants_path : str
        Path to the CSV file containing ammonia plant data.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the ammonia plant data with columns: location, capacity, x, y, country, etc.
    """
    # Load plants from provided CSV
    ammonia_plants = read_csv_nafix(ammonia_plants_path)
    ammonia_plants = ammonia_plants.rename(
        columns={
            "Latitude": "y",
            "Longitude": "x",
            "Plant": "plant",
            "Ammonia [kt/a]": "capacity",
        }
    )

    # Convert country names to ISO2 codes
    cc = coco.CountryConverter()
    ammonia_plants["country"] = cc.convert(ammonia_plants["Country"], to="ISO2")

    # Select relevant columns and set index
    ammonia_plants = ammonia_plants[["plant", "capacity", "x", "y", "country"]]
    ammonia_plants.index.name = "ktonNH3/a"

    return ammonia_plants


def combine_ammonia_plants(
    plants: pd.DataFrame, us_plants: pd.DataFrame
) -> pd.DataFrame:
    """
    Combine custom and US ammonia plant data into a single DataFrame.

    Parameters
    ----------
    plants : pd.DataFrame
        DataFrame containing custom ammonia plant data.
    us_plants : pd.DataFrame
        DataFrame containing US ammonia plant data.

    Returns
    -------
    pd.DataFrame
        A combined DataFrame containing ammonia plant data for both US and custom plants.
    """
    # Combine US and custom plant data
    ammonia_plants = pd.concat([plants, us_plants], ignore_index=True)

    # Drop rows with missing capacity
    ammonia_plants = ammonia_plants.dropna(subset=["capacity"])

    return ammonia_plants


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_ammonia_production")

    configure_logging(snakemake)

    # Load inputs
    ammonia_plants_path = snakemake.input.ammonia_plants
    us_cities_path = snakemake.input.us_cities
    usgs_ammonia_dataset_path = snakemake.input.usgs_ammonia_dataset

    # Extract country-level production totals
    ammonia_totals = extract_country_level_production(usgs_ammonia_dataset_path)

    # Extract ammonia plant data
    ammonia_plants = extract_plant_data(ammonia_plants_path)

    # Load US cities data
    us_cities = read_csv_nafix(us_cities_path)

    # Extract US plant data
    us_ammonia_plants = extract_us_plant_data(usgs_ammonia_dataset_path, us_cities)

    # Combine ammonia plants data
    combined_ammonia_plants = combine_ammonia_plants(ammonia_plants, us_ammonia_plants)

    # Save the results
    to_csv_nafix(ammonia_totals, snakemake.output.ammonia_production)
    to_csv_nafix(combined_ammonia_plants, snakemake.output.ammonia_plants, index=False)
