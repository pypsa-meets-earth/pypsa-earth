# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Process EDGAR CO2 emission data from the raw Excel file into a clean CSV.

Reads the EDGAR CO2 fossil fuel emission Excel file, filters to
'Public electricity and heat production' entries, and saves the result
as a CSV file indexed by ISO3 country code for use by prepare_network.
"""

import pandas as pd
from _helpers import configure_logging, create_logger

logger = create_logger(__name__)


def process_edgar_emission_data(excel_file):
    """
    Process EDGAR CO2 emission Excel data into a clean DataFrame.

    Detects the CO2 fossil fuel emissions sheet automatically, filters to
    'Public electricity and heat production' rows, and returns a DataFrame
    indexed by Country_code_A3 with year columns (Y_YYYY).

    Parameters
    ----------
    excel_file : str
        Path to the EDGAR CO2 emission Excel file (.xls or .xlsx).

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by Country_code_A3 with columns Y_YYYY for each
        available year.
    """
    xl = pd.ExcelFile(excel_file)
    sheet_names = xl.sheet_names

    target_sheet = None
    for sheet in sheet_names:
        if "EM_CO2_fossil" in sheet:
            target_sheet = sheet
            break

    if target_sheet is None:
        raise ValueError(
            "Could not find CO2 fossil emissions sheet in EDGAR data. "
            f"Available sheets: {sheet_names}"
        )

    logger.info(f"Processing EDGAR emission data from sheet: '{target_sheet}'")

    df = pd.read_excel(excel_file, sheet_name=target_sheet, skiprows=8)
    df.columns = df.iloc[0]
    df = df.iloc[1:].set_index("Country_code_A3")
    df = df.loc[
        df["IPCC_for_std_report_desc"] == "Public electricity and heat production"
    ]

    year_cols = [
        col for col in df.columns if isinstance(col, str) and col.startswith("Y_")
    ]
    df = df[year_cols].astype(float).ffill(axis=1).bfill(axis=1)

    logger.info(
        f"Processed emission data for {len(df)} countries, "
        f"years {year_cols[0]}–{year_cols[-1]}."
    )
    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_co2_emissions")

    configure_logging(snakemake)

    df = process_edgar_emission_data(snakemake.input.edgar)
    df.to_csv(snakemake.output.emissions)
    logger.info(f"Emission data saved to '{snakemake.output.emissions}'.")
