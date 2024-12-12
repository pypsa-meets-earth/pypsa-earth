# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import math

import country_converter as coco
import numpy as np
import pandas as pd
import pycountry
import requests
from _helpers import content_retrieve
from geopy.geocoders import Nominatim


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
                    print(
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


def create_steel_db():
    # Global Steel Plant Tracker data set you requested from Global Energy Monitor from the link below:

    # The following excel file was downloaded from the following webpage
    # https://globalenergymonitor.org/wp-content/uploads/2023/03/Global-Steel-Plant-Tracker-2023-03.xlsx . The dataset contains 1433 Steel plants globally.

    url = "https://globalenergymonitor.org/wp-content/uploads/2023/03/Global-Steel-Plant-Tracker-2023-03.xlsx"

    df_steel = pd.read_excel(
        content_retrieve(url),
        index_col=0,
        sheet_name="Steel Plants",
        header=0,
    )

    df_steel = df_steel[
        [
            "Plant name (English)",
            "Country",
            "Coordinates",
            "Coordinate accuracy",
            "Status",
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
    df_steel = df_steel.loc[df_steel["Status"] == "operating"]

    # Create a column with iso2 country code
    cc = coco.CountryConverter()
    Country = pd.Series(df_steel["Country"])
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
    df_steel["Main production process"].value_counts()

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
    df_steel["Main production process"].value_counts()

    # Dict to map the technology names of the source to that expected in the workflow
    iron_techs = {
        "electric": "Electric arc",
        "integrated (BF)": "Integrated steelworks",
        "integrated (DRI)": "DRI + Electric arc",
        "ironmaking (BF)": "Integrated steelworks",
        "ironmaking (DRI)": "DRI + Electric arc",
        "oxygen": "Integrated steelworks",
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
    df_steel.capacity = pd.to_numeric(df_steel.capacity)
    df_steel["technology"] = df_steel["Main production process"].apply(
        lambda x: iron_techs[x]
    )
    df_steel.x = df_steel.x.apply(lambda x: eval(x))
    df_steel.y = df_steel.y.apply(lambda y: eval(y))

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


def create_cement_db():
    # -------------
    # CEMENT
    # -------------
    # The following excel file was downloaded from the following webpage https://www.cgfi.ac.uk/spatial-finance-initiative/geoasset-project/cement/.
    # The dataset contains 3117 cement plants globally.
    fn = "https://www.cgfi.ac.uk/wp-content/uploads/2021/08/SFI-Global-Cement-Database-July-2021.xlsx"
    storage_options = {"User-Agent": "Mozilla/5.0"}
    cement_orig = pd.read_excel(
        fn,
        index_col=0,
        storage_options=storage_options,
        sheet_name="SFI_ALD_Cement_Database",
        header=0,
    )

    df_cement = cement_orig.copy()
    df_cement = df_cement[
        [
            "country",
            "iso3",
            "latitude",
            "longitude",
            "status",
            "plant_type",
            "capacity",
            "year",
            "city",
        ]
    ]
    df_cement = df_cement.rename(
        columns={
            "country": "Country",
            "latitude": "y",
            "longitude": "x",
            "city": "location",
        }
    )
    df_cement["unit"] = "Kt/yr"
    df_cement["technology"] = "Cement"
    df_cement["capacity"] = df_cement["capacity"] * 1000
    # Keep only operating steel plants
    df_cement = df_cement.loc[df_cement["status"] == "Operating"]

    # Create a column with iso2 country code
    cc = coco.CountryConverter()
    iso3 = pd.Series(df_cement["iso3"])
    df_cement["country"] = cc.pandas_convert(series=iso3, to="ISO2")

    # Dropping the null capacities reduces the dataframe from 3000+  rows to 1672 rows
    na_index = df_cement[df_cement.capacity.isna()].index
    print(
        "There are {} out of {} total cement plants with unknown capacities, setting value to country average".format(
            len(na_index), len(df_cement)
        )
    )
    avg_c_cap = df_cement.groupby(df_cement.country)["capacity"].mean()
    df_cement["capacity"] = df_cement.apply(
        lambda x: (
            avg_c_cap[x["country"]] if math.isnan(x["capacity"]) else x["capacity"]
        ),
        axis=1,
    )

    df_cement["quality"] = "actual"
    df_cement.loc[na_index, "quality"] = "actual"  # TODO change

    df_cement = df_cement.reset_index()
    df_cement = df_cement.rename(columns={"uid": "ID"})
    df_cement.capacity = pd.to_numeric(df_cement.capacity)

    return df_cement[
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


def create_refineries_df():
    # -------------
    # OIL REFINERIES
    # -------------
    # The data were downloaded directly from arcgis server using a query found on this webpage:
    # https://www.arcgis.com/home/item.html?id=a6979b6bccbf4e719de3f703ea799259&sublayer=0#data
    # and https://www.arcgis.com/home/item.html?id=a917ac2766bc47e1877071f0201b6280

    # The dataset contains 536 global Oil refineries.

    base_url = "https://services.arcgis.com"
    facts = "/jDGuO8tYggdCCnUJ/arcgis/rest/services/Global_Oil_Refinery_Complex_and_Daily_Capacity/FeatureServer/0/query?f=json&where=1%3D1&returnGeometry=false&spatialRel=esriSpatialRelIntersects&outFields=*&orderByFields=FID%20ASC&resultOffset=0&resultRecordCount=537&cacheHint=true&quantizationParameters=%7B%22mode%22%3A%22edit%22%7D"

    first_response = requests.get(base_url + facts)
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


def create_paper_df():
    # -------------
    # Paper
    # -------------
    # The following excel file was downloaded from the following webpage https://www.cgfi.ac.uk/spatial-finance-initiative/geoasset-project/cement/ . The dataset contains 3117 cement plants globally.

    fn = "https://www.cgfi.ac.uk/wp-content/uploads/2023/03/SFI_ALD_Pulp_Paper_Sample_LatAm_Jan_2023.xlsx"

    storage_options = {"User-Agent": "Mozilla/5.0"}
    paper_orig = pd.read_excel(
        fn,
        index_col=0,
        storage_options=storage_options,
        sheet_name="SFI_ALD_PPM_LatAm",
        header=0,
    )

    df_paper = paper_orig.copy()
    df_paper = df_paper[
        [
            "country",
            "iso3",
            "latitude",
            "longitude",
            "status",
            "primary_product",
            "capacity_paper",
            "city",
        ]
    ]

    df_paper = df_paper.rename(
        columns={
            "country": "Country",
            "latitude": "y",
            "longitude": "x",
            "city": "location",
            "capacity_paper": "capacity",
        }
    )
    df_paper["unit"] = "10kt/yr"
    df_paper["technology"] = "Paper"
    df_paper["capacity"] = df_paper["capacity"]

    df_paper.capacity = df_paper.capacity.apply(
        lambda x: x if type(x) == int or type(x) == int == float else np.nan
    )

    # Keep only operating steel plants
    # df_paper = df_paper.loc[df_paper["status"] == "Operating"]

    # Create a column with iso2 country code
    cc = coco.CountryConverter()
    iso3 = pd.Series(df_paper["iso3"])
    df_paper["country"] = cc.pandas_convert(series=iso3, to="ISO2")

    # Dropping the null capacities reduces the dataframe from 3000+  rows to 1672 rows
    na_index = df_paper[df_paper.capacity.isna()].index
    print(
        "There are {} out of {} total paper plants with unknown capacities, setting value to country average".format(
            len(na_index), len(df_paper)
        )
    )
    avg_c_cap = df_paper.groupby(df_paper.country)["capacity"].mean()
    na_index

    df_paper["capacity"] = df_paper.apply(
        lambda x: (
            avg_c_cap[x["country"]] if math.isnan(x["capacity"]) else x["capacity"]
        ),
        axis=1,
    )

    df_paper["quality"] = "actual"
    df_paper.loc[na_index, "quality"] = "actual"  # TODO change
    df_paper.capacity = pd.to_numeric(df_paper.capacity)

    df_paper = df_paper.reset_index()
    df_paper = df_paper.rename(columns={"uid": "ID"})

    industrial_database_paper = df_paper[
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

    no_infp_index = industrial_database_paper[
        industrial_database_paper.y == "No information"
    ].index
    print(
        "Setting plants of countries with no values for paper plants to 1.0".format(
            len(na_index), len(df_paper)
        )
    )
    industrial_database_paper = industrial_database_paper.drop(no_infp_index)
    industrial_database_paper.capacity = industrial_database_paper.capacity.fillna(1)

    return industrial_database_paper


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industrial_database",
            simpl="",
            clusters="4",
            ll="c1",
            opts="Co2L-4H",
            planning_horizons="2030",
            sopts="144H",
            discountrate=0.071,
            demand="AB",
        )

    industrial_database_steel = create_steel_db()
    industrial_database_cement = create_cement_db()
    industrial_database_refineries = create_refineries_df()
    industrial_database_paper = create_paper_df()

    industrial_database = pd.concat(
        [
            industrial_database_steel,
            industrial_database_cement,
            industrial_database_refineries,
            industrial_database_paper,
        ]
    )

    industrial_database.to_csv(
        snakemake.output["industrial_database"], header=True, index=0
    )
