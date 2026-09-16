# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import math

import country_converter as coco
import numpy as np
import pandas as pd
from _helpers import read_csv_nafix


def create_cement_db(fn):
    """
    Read SFI cement database
    """
    cement_orig = pd.read_excel(
        fn,
        index_col=0,
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


def create_paper_df(fn):
    """
    Pre-process cement database
    """

    paper_orig = pd.read_excel(
        fn,
        index_col=0,
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

    cc = coco.CountryConverter()
    iso3 = pd.Series(df_paper["iso3"])
    df_paper["country"] = cc.pandas_convert(series=iso3, to="ISO2")

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


def create_ammonia_db(ammonia_plants_file: str) -> pd.DataFrame:
    """
    Read ammonia plants database from resources.

    The ammonia_plants.csv file is created by build_ammonia_production.py
    and contains combined US and EU plant data with coordinates.

    Parameters
    ----------
    ammonia_plants_file : str
        Path to the ammonia plants CSV file.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing ammonia plant information with columns:
        ['country', 'y', 'x', 'location', 'technology', 'capacity', 'unit', 'quality', 'ID'].
    """
    # Load ammonia plants data
    df_ammonia = read_csv_nafix(ammonia_plants_file)

    # Set location to plant name
    df_ammonia["location"] = df_ammonia["plant"]

    # Set technology to Haber-Bosch (the standard ammonia synthesis process)
    df_ammonia["technology"] = "Haber-Bosch"

    # Unit is kt/yr (kilotons per annum) of NH3
    df_ammonia["unit"] = "kt/yr"

    # Quality is exact (from actual plant data)
    df_ammonia["quality"] = "exact"

    # Use plant index as ID
    df_ammonia["ID"] = df_ammonia.index

    return df_ammonia[
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

        snakemake = mock_snakemake(
            "build_industrial_database",
            simpl="",
            clusters="4",
            ll="c1",
            opts="Co2L-4H",
            planning_horizons="2030",
            sopts="144H",
            discountrate=0.071,
        )

    # Load parameters
    ammonia_plants_file = snakemake.input.ammonia_plants

    industrial_database_steel = read_csv_nafix(snakemake.input.steel_raw)
    industrial_database_cement = create_cement_db(snakemake.params.url_cement)
    industrial_database_refineries = read_csv_nafix(snakemake.input.refineries_raw)
    industrial_database_paper = create_paper_df(snakemake.params.url_paper)
    industrial_database_ammonia = create_ammonia_db(ammonia_plants_file)

    industrial_database = pd.concat(
        [
            industrial_database_steel,
            industrial_database_cement,
            industrial_database_refineries,
            industrial_database_paper,
            industrial_database_ammonia,
        ]
    )

    industrial_database.to_csv(
        snakemake.output["industrial_database"], header=True, index=0
    )
