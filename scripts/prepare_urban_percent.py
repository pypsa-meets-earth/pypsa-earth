# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
import country_converter as coco
import pandas as pd
from _helpers import read_csv_nafix

# from _helpers import configure_logging


# logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_urban_percent")

    df = read_csv_nafix(snakemake.input.urban_percent_raw).copy()

    # Select the columns that we need to keep
    df = df[
        [
            "Year",
            "Economy Label",
            "Urban population as percentage of total population",
        ]
    ]

    # Keep only years above 2020
    df = df.loc[(df["Year"] >= 2020)]

    # Add ISO2 country code for each country
    cc = coco.CountryConverter()
    Economy_Label = pd.Series(df["Economy Label"])
    df["country"] = cc.pandas_convert(
        series=Economy_Label, to="ISO2", not_found="not found"
    )

    # Drop isos that were not found:
    df = df.loc[df["country"] != "not found"]

    # Drop region names where country column contains list of countries
    df = df.loc[df.country.apply(lambda x: isinstance(x, str)), :]

    # Reduce the data to one value for the urban percent per country and year
    df = df.groupby(["country", "Year"], as_index=False).mean(numeric_only=True)

    df = df.set_index("country")

    # Save
    df.to_csv(snakemake.output[0], sep=",", encoding="utf-8", header="true")
