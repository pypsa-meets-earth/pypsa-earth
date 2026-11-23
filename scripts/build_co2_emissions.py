# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Process global EDGAR CO2 emission data from the raw Excel file into a clean CSV.

Reads the EDGAR CO2 fossil fuel emission Excel file, filters to
'Public electricity and heat production' entries, and saves the result
as a CSV with country identifier columns (country_code_a3, country_code_a2,
country_name) and year columns for use by prepare_network.
"""

import pandas as pd
from _helpers import configure_logging, create_logger, three_2_two_digits_countries

logger = create_logger(__name__)


def process_edgar_emission_data(
    excel_file: str,
    target_sheet: str = "v6.0_EM_CO2_fossil_IPCC2006",
) -> pd.DataFrame:
    """
    Process EDGAR CO2 emission Excel data into a clean DataFrame.

    Detects the CO2 fossil fuel emissions sheet automatically, filters to
    'Public electricity and heat production' rows, and returns a DataFrame
    with country identifier columns (country_code_a3, country_code_a2, country_name)
    and year columns (Y_YYYY).

    Parameters
    ----------
    excel_file : str
        Path to the EDGAR CO2 emission Excel file (.xls or .xlsx).
    target_sheet : str, optional
        Name of the sheet in the Excel file to process. Defaults to
        'v6.0_EM_CO2_fossil_IPCC2006'.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns country_code_a3 (three-letter ISO code),
        country_code_a2 (two-letter ISO code), country_name (full country name),
        and Y_YYYY columns for each available year.
    """
    logger.info(f"Processing EDGAR emission data from sheet: '{target_sheet}'")

    df = pd.read_excel(excel_file, sheet_name=target_sheet, skiprows=8)
    df.columns = df.iloc[0]
    df = df.iloc[1:]
    df = df.loc[
        df["ipcc_code_2006_for_standard_report_name"]
        == "Main Activity Electricity and Heat Production"
    ]  # TODO: extend filter to other relevant sectors

    year_cols = [
        col for col in df.columns if isinstance(col, str) and col.startswith("Y_")
    ]
    df = df[["Country_code_A3", "Name"] + year_cols].copy()
    df = df.rename(
        columns={
            "Country_code_A3": "country_code_a3",
            "Name": "country_name",
        }
    )
    df["country_code_a2"] = three_2_two_digits_countries(df["country_code_a3"])
    df = df[["country_code_a3", "country_code_a2", "country_name"] + year_cols]
    df[year_cols] = df[year_cols].astype(float).ffill(axis=1).bfill(axis=1)

    # rename year columns to remove 'Y_' prefix
    df = df.rename(columns={col: col[2:] for col in year_cols})

    logger.info(
        f"Processed emission data for {len(df)} countries, "
        f"years {year_cols[0][2:]}–{year_cols[-1][2:]}."
    )
    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_co2_emissions")

    configure_logging(snakemake)

    df = process_edgar_emission_data(snakemake.input.edgar)
    df.to_csv(snakemake.output.emissions, index=False)
    logger.info(f"Emission data saved to '{snakemake.output.emissions}'.")
