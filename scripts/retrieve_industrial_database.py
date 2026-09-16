# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Retrieve the steel plant and oil refinery datasets used by
build_industrial_database.py.
"""

import country_converter as coco
import pandas as pd
import pycountry
import requests
from _helpers import content_retrieve
from geopy.geocoders import Nominatim

logger = create_logger(__name__)


def get_cocode_from_name(df, country_column_name):
    country_codes = {}

    for country in pycountry.countries:
        country_codes[country.name] = country.alpha_2

    df["country"] = df[country_column_name].map(country_codes)
    return df


def get_cocode_from_coords(df):
    geolocator = Nominatim(user_agent="geoapi")  # Initialize geolocator

    # Initialize an empty list to store country codes
    country_codes = []

    for index, row in df.iterrows():
        # Get latitude and longitude from the row
        latitude = row["Latitude"]
        longitude = row["Longitude"]

        # Perform reverse geocoding to get location information
        tries = 0
        location = None
        while tries < 10:
            try:
                location = geolocator.reverse((latitude, longitude), exactly_one=True)
                break
            except:
                tries += 1
                if tries == 10:
                    logger.error(
                        "Country code of location ({},{}) could not be geocoded after 10 tries.".format(
                            latitude, longitude
                        )
                    )

        if location and location.raw.get("address", {}).get("country_code"):
            # Extract and append the country code to the list
            country_code = location.raw["address"]["country_code"].upper()
            country_codes.append(country_code)
        else:
            country_codes.append(None)

    # Add the country code list as a new column to the DataFrame
    df["country"] = country_codes

    return df


def create_steel_db(fn):
    df_steel = pd.read_excel(
        content_retrieve(fn),
        index_col=0,
        sheet_name="Steel Plants",
        header=0,
    )

    df_steel = df_steel[
        [
            "Plant name (English)",
            "Country/Area",
            "Coordinates",
            "Coordinate accuracy",
            "Capacity operating status",
            "Start date",
            "Plant age (years)",
            "Nominal crude steel capacity (ttpa)",
            "Nominal BOF steel capacity (ttpa)",
            "Nominal EAF steel capacity (ttpa)",
            "Nominal OHF steel capacity (ttpa)",
            "Nominal iron capacity (ttpa)",
            "Nominal BF capacity (ttpa)",
            "Nominal DRI capacity (ttpa)",
            "Ferronickel capacity (ttpa)",
            "Sinter plant capacity (ttpa)",
            "Coking plant capacity (ttpa)",
            "Pelletizing plant capacity (ttpa)",
            "Category steel product",
            "Main production process",
            "Municipality",
        ]
    ]

    # Keep only operating steel plants
    df_steel = df_steel.loc[df_steel["Capacity operating status"] == "operating"]

    # Create a column with iso2 country code
    cc = coco.CountryConverter()
    Country = pd.Series(df_steel["Country/Area"])
    df_steel["country"] = cc.pandas_convert(series=Country, to="ISO2")

    # Split Coordeinates column into x and y columns
    df_steel[["y", "x"]] = df_steel["Coordinates"].str.split(",", expand=True)

    # Drop Coordinates column as it contains a ',' and is not needed anymore
    df_steel = df_steel.drop(columns="Coordinates", axis=1)

    # Fetch steel plants that uses DRI and BF techs and drop them from main df
    mixed_steel_plants = df_steel[
        df_steel["Main production process"] == "integrated (BF and DRI)"
    ].copy()
    df_steel = df_steel.drop(mixed_steel_plants.index)

    # Separate the two techs in two dataframes
    DRI_share = mixed_steel_plants.copy()
    BF_share = mixed_steel_plants.copy()
    BF_share["Main production process"] = "integrated (BF)"
    DRI_share["Main production process"] = "integrated (DRI)"

    # Calculate the share of both techs according to the capacities of iron production
    BF_share["Nominal crude steel capacity (ttpa)"] = BF_share[
        "Nominal crude steel capacity (ttpa)"
    ] * mixed_steel_plants.apply(
        lambda x: x["Nominal BF capacity (ttpa)"] / x["Nominal iron capacity (ttpa)"],
        axis=1,
    )
    DRI_share["Nominal crude steel capacity (ttpa)"] = (
        mixed_steel_plants["Nominal crude steel capacity (ttpa)"]
        - BF_share["Nominal crude steel capacity (ttpa)"]
    )

    # Add suffix to the index to differentiate between them in the main df
    DRI_share.index += "_DRI"
    BF_share.index += "_BF"

    # Merge them back to the main df
    df_steel = pd.concat([df_steel, BF_share, DRI_share])

    # Remove plants with unknown production technology
    unknown_ind = df_steel[
        df_steel["Main production process"].str.contains("unknown")
    ].index
    df_steel = df_steel.drop(unknown_ind)
    if len(unknown_ind) > 0:
        print(
            "dropped {0} steel/iron plants with unknown production technology of total {1} plants".format(
                len(unknown_ind), len(df_steel)
            )
        )

    # Dict to map the technology names of the source to that expected in the workflow
    iron_techs = {
        "electric": "Electric arc",
        "integrated (BF)": "Integrated steelworks",
        "integrated (DRI)": "DRI + Electric arc",
        "ironmaking (BF)": "Integrated steelworks",
        "ironmaking (DRI)": "DRI + Electric arc",
        "oxygen": "Integrated steelworks",
        "ironmaking (other)": "Integrated steelworks",
        "steelmaking (other)": "Integrated steelworks",
        "electric, oxygen": "Electric arc",
    }

    # Creating the necessary columns in the dataframe
    iron_making = df_steel[
        df_steel["Main production process"].str.contains("ironmaking")
    ].index
    df_steel.loc[iron_making, "Nominal crude steel capacity (ttpa)"] = df_steel.loc[
        iron_making, "Nominal iron capacity (ttpa)"
    ]
    df_steel["unit"] = "kt/yr"
    df_steel["quality"] = "exact"
    df_steel = df_steel.reset_index()
    df_steel = df_steel.rename(
        columns={
            "Nominal crude steel capacity (ttpa)": "capacity",
            "Municipality": "location",
            "Plant ID": "ID",
        }
    )
    df_steel["technology"] = df_steel["Main production process"].replace(iron_techs)

    for col in ["capacity", "x", "y"]:
        df_steel[col] = pd.to_numeric(df_steel[col], errors="coerce")

    return df_steel[
        [
            "country",
            "y",
            "x",
            "location",
            "technology",
            "capacity",
            "unit",
            "quality",
            "ID",
        ]
    ].dropna()


def create_refineries_df(fn):
    """
    Pre-process refineries data
    """
    first_response = requests.get(fn)
    response_list = first_response.json()

    data = []
    for response in response_list["features"]:
        data.append(
            {
                "FID_": response["attributes"].get("FID_"),
                "Company": response["attributes"].get("Company"),
                "Name": response["attributes"].get("Name"),
                "City": response["attributes"].get("City"),
                "Facility": response["attributes"].get("Facility"),
                "Prov_State": response["attributes"].get("Prov_State"),
                "Country": response["attributes"].get("Country"),
                "Address": response["attributes"].get("Address"),
                "Zip": response["attributes"].get("Zip"),
                "County": response["attributes"].get("County"),
                "PADD": response["attributes"].get("PADD"),
                "Capacity": response["attributes"].get("Capacity"),
                "Longitude": response["attributes"].get("Longitude"),
                "Latitude": response["attributes"].get("Latitude"),
                "Markets": response["attributes"].get("Markets"),
                "CORPORATIO": response["attributes"].get("CORPORATIO"),
            }
        )

    df = pd.DataFrame(data)

    df = get_cocode_from_name(df, "Country")

    df_nans = df[df.country.isna()]
    df = df.dropna(axis=0)

    df_bylocation = get_cocode_from_coords(df_nans)

    df_refineries = pd.concat([df, df_bylocation])

    # Creating the necessary columns in the dataframe
    # df_refineries["technology"] = df_refineries["Main production process"].apply(lambda x: iron_techs[x])
    df_refineries["unit"] = "bpd"
    df_refineries["quality"] = "exact"
    df_refineries["technology"] = "HVC"

    df_refineries = df_refineries.rename(
        columns={
            "Capacity": "capacity",
            "Prov_State": "location",
            "Latitude": "y",
            "Longitude": "x",
            "FID_": "ID",
        }
    )
    df_refineries = df_refineries.reset_index()
    df_refineries.capacity = pd.to_numeric(df_refineries.capacity)

    return df_refineries[
        [
            "country",
            "y",
            "x",
            "location",
            "technology",
            "capacity",
            "unit",
            "quality",
            "ID",
        ]
    ]


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_industrial_database")

    industrial_database_steel = create_steel_db(snakemake.params.url_steel)
    industrial_database_steel.to_csv(snakemake.output.steel_raw, index=False)

    industrial_database_refineries = create_refineries_df(
        snakemake.params.url_refineries
    )
    industrial_database_refineries.to_csv(snakemake.output.refineries_raw, index=False)
