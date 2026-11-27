# SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Shared UNSD SDMX data loading and normalization."""

from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd
from download_unsd_sdmx import country_to_unsd_area


def read_sdmx_labels(root, codelist_id):
    """Map official SDMX codes to their English descriptions."""
    ns = {
        "s": "http://www.sdmx.org/resources/sdmxml/schemas/v2_1/structure",
        "c": "http://www.sdmx.org/resources/sdmxml/schemas/v2_1/common",
    }

    codelist = root.find(
        f".//s:Codelist[@id='{codelist_id}']",
        ns,
    )

    if codelist is None:
        raise ValueError(f"Missing UNSD codelist: {codelist_id}")

    return {
        code.get("id"): code.findtext("c:Name", namespaces=ns)
        for code in codelist.findall("s:Code", ns)
    }


def load_unsd(unsd_dir, countries):
    """Read cached SDMX observations using the official UNSD classifications."""
    root = ET.parse(Path(unsd_dir) / "structure.xml").getroot()

    commodities = read_sdmx_labels(root, "CL_COMMODITY_NRG")
    transactions = read_sdmx_labels(root, "CL_TRANSACTION_NRG")

    frames = []

    for country in countries:
        code = country_to_unsd_area(country)
        source = Path(unsd_dir) / f"{code}.csv"

        df = pd.read_csv(
            source,
            dtype={
                "COMMODITY": str,
                "TRANSACTION": str,
            },
            low_memory=False,
        )

        df = df[df["FREQ"].eq("A")].copy()

        df["country"] = country
        df["Commodity"] = df["COMMODITY"].map(commodities)

        # Normalize SDMX labels to the existing PyPSA-Earth fuel mapping.
        commodity_aliases = {
            "7000": "Electricity",
            "4680": "Fuel oil",
            "0100": "Hard coal",
            "0200": "Brown coal",
            "4630": "Liquefied petroleum gas (LPG)",
            "5140": "Black Liquor",
            "5150": "Vegetal waste",
            "0379": "Recovered gases",
        }

        df["Commodity"] = df["COMMODITY"].map(commodity_aliases).fillna(df["Commodity"])
        df["Commodity"] = df["Commodity"].replace(
            {"Aviation Gasoline": "Aviation gasoline"}
        )
        df["Transaction"] = df["TRANSACTION"].map(transactions)
        df["transaction_code"] = df["TRANSACTION"]
        df["Year"] = pd.to_numeric(df["TIME_PERIOD"], errors="coerce")

        # Apply the SDMX unit multiplier before unit conversion.
        multiplier = pd.to_numeric(df["UNIT_MULT"], errors="coerce")
        value = pd.to_numeric(df["OBS_VALUE"], errors="coerce")

        df["Quantity"] = value * (10.0**multiplier)

        unit = df["UNIT_MEASURE"].astype(str).str.upper()

        units = {
            "TN": "Metric tons,  thousand",
            "M3": "Cubic metres, thousand",
            "TJ": "Terajoules",
            "GWHR": "Kilowatt-hours, million",
            "KWH": "Kilowatt-hours, million",
            "MWH": "Kilowatt-hours, million",
            "GWH": "Kilowatt-hours, million",
        }

        df["Unit"] = unit.map(units)

        # PyPSA-Earth conversion factors are expressed per thousand tonnes
        # or thousand cubic metres.
        mass_or_volume = unit.isin(["TN", "M3"])
        df.loc[mass_or_volume, "Quantity"] /= 1000

        # Normalize electricity observations to million kWh.
        df.loc[unit.eq("KWH"), "Quantity"] /= 1e6
        df.loc[unit.eq("MWH"), "Quantity"] /= 1e3

        df["Commodity - Transaction"] = (
            df["Commodity"].astype(str) + " - " + df["Transaction"].astype(str)
        )

        df = df.dropna(subset=["Commodity", "Transaction", "Year"])

        df["Year"] = df["Year"].astype(int)

        frames.append(df)

    return pd.concat(frames, ignore_index=True)


SECTOR_TRANSACTIONS = {
    "consumption by households": ("1231",),
    "services": ("1235",),
    "agriculture": ("1232",),
    "road": ("1221",),
    "rail": ("1222",),
    "aviation": ("1223", "052"),
    "navigation": ("1224", "051"),
    "other energy": ("1234",),
    "non energy use": ("11",),
}


def select_base_year(df, year, max_year_gap=0):
    """Select observations without using data from future years.

    Return selected observations and historical series excluded
    because their latest observation exceeds the allowed year gap.
    """
    if max_year_gap < 0:
        raise ValueError("max_year_gap must be non-negative")

    keys = ["country", "COMMODITY", "transaction_code"]

    eligible = df.loc[df.Year.le(year) & df.Quantity.notna()].copy()

    latest_year = eligible.groupby(keys)["Year"].transform("max")
    latest = eligible.loc[eligible.Year.eq(latest_year)].copy()

    latest["year_gap"] = year - latest.Year

    selected = latest.loc[latest.year_gap.le(max_year_gap)].copy()
    missing = latest.loc[latest.year_gap.gt(max_year_gap)].copy()

    if selected.duplicated(keys).any():
        raise ValueError("Multiple UNSD observations for the same series")

    return selected, missing


def select_sector(df, sector):
    """Select final consumption using official UNSD transaction codes."""
    return df.loc[df.transaction_code.isin(SECTOR_TRANSACTIONS[sector])].copy()


INDUSTRY_TRANSACTIONS = {
    "1211": "iron and steel",
    "1213": "chemical and petrochemical",
    "1214A": "non-ferrous metals",
    "1214B": "non-metallic minerals",
    "1214C": "transport equipment",
    "1214D": "machinery",
    "1214E": "mining and quarrying",
    "1214F": "food and tobacco",
    "1214G": "paper pulp and print",
    "1214H": "wood and wood products",
    "1214I": "construction",
    "1214J": "textile and leather",
    "1214O": "other",
}


def select_industry(df):
    """Select disaggregated industrial consumption without sector totals."""
    rows = df.loc[df.transaction_code.isin(INDUSTRY_TRANSACTIONS)].copy()

    rows["Transaction"] = rows.transaction_code.map(INDUSTRY_TRANSACTIONS)

    return rows


INDUSTRY_TRANSACTIONS = {
    "1211": "iron and steel",
    "1213": "chemical and petrochemical",
    "1214A": "non-ferrous metals",
    "1214B": "non-metallic minerals",
    "1214C": "transport equipment",
    "1214D": "machinery",
    "1214E": "mining and quarrying",
    "1214F": "food and tobacco",
    "1214G": "paper pulp and print",
    "1214H": "wood and wood products",
    "1214I": "construction",
    "1214J": "textile and leather",
    "1214O": "other",
}


def select_industry(df):
    """Select disaggregated industrial consumption without sector totals."""
    rows = df.loc[df.transaction_code.isin(INDUSTRY_TRANSACTIONS)].copy()

    rows["Transaction"] = rows.transaction_code.map(INDUSTRY_TRANSACTIONS)

    return rows


def read_commodity_parents(structure_path):
    """Read the official UNSD commodity hierarchy."""
    root = ET.parse(structure_path).getroot()
    ns = {
        "s": "http://www.sdmx.org/resources/sdmxml/schemas/v2_1/structure",
    }

    codes = root.findall(".//s:Codelist[@id='CL_COMMODITY_NRG']/s:Code", ns)

    if not codes:
        raise ValueError("UNSD commodity codelist not found")

    parents = {}

    for code in codes:
        parent = code.find("./s:Parent/Ref", ns)
        if parent is not None:
            parents[code.get("id")] = parent.get("id")

    return parents


def reconcile_commodities(rows, parents):
    """Remove descendants already represented by aggregate commodities."""
    rows = rows.reset_index(drop=True)
    discarded = set()

    for _, group in rows.groupby(["country", "transaction_code"]):
        available = set(group.COMMODITY)

        for index, code in group.COMMODITY.items():
            ancestor = parents.get(code)
            visited = set()

            while ancestor and ancestor not in visited:
                visited.add(ancestor)

                if ancestor in available:
                    discarded.add(index)
                    break

                ancestor = parents.get(ancestor)

    return rows.drop(index=list(discarded)).copy()


def convert_to_twh(rows, factors):
    """Convert normalized UNSD observations to TWh."""
    rows = rows.copy()

    mass_or_volume = rows.Unit.isin(
        ["Metric tons,  thousand", "Cubic metres, thousand"]
    )
    electricity = rows.Unit.eq("Kilowatt-hours, million")
    energy = rows.Unit.eq("Terajoules")

    rows["Quantity_TWh"] = float("nan")

    rows.loc[mass_or_volume, "Quantity_TWh"] = rows.loc[
        mass_or_volume, "Quantity"
    ] * rows.loc[mass_or_volume, "Commodity"].map(factors)

    rows.loc[electricity, "Quantity_TWh"] = rows.loc[electricity, "Quantity"] / 1000

    rows.loc[energy, "Quantity_TWh"] = rows.loc[energy, "Quantity"] / 3600

    return rows


def prepare_demand_rows(
    directory, countries, year, factors, commodities, max_year_gap=0, industry=False
):
    """Prepare normalized and reconciled UNSD demand observations."""
    import logging

    logger = logging.getLogger(__name__)
    rows = load_unsd(directory, countries)

    if industry:
        codes = set(INDUSTRY_TRANSACTIONS)
    else:
        codes = {
            code
            for transactions in SECTOR_TRANSACTIONS.values()
            for code in transactions
        }

    rows = rows.loc[rows.transaction_code.isin(codes)]

    selected, missing = select_base_year(rows, year, max_year_gap)

    if not missing.empty:
        logger.warning(
            "%s historical UNSD series excluded by year-gap policy.",
            len(missing),
        )

    unmapped = selected.loc[~selected.Commodity.isin(commodities)]

    if not unmapped.empty:
        logger.warning(
            "Unmapped UNSD commodities: %s",
            sorted(unmapped.Commodity.unique()),
        )

    selected = selected.loc[selected.Commodity.isin(commodities)]
    selected = convert_to_twh(selected, factors)

    invalid = selected.Quantity_TWh.isna()
    if invalid.any():
        raise ValueError(
            "Missing UNSD energy conversions: "
            + str(selected.loc[invalid, "Commodity"].unique())
        )

    parents = read_commodity_parents(Path(directory) / "structure.xml")

    selected = reconcile_commodities(selected, parents)

    missing_countries = set(countries) - set(selected.country)

    if missing_countries:
        raise ValueError(f"No UNSD demand observations for {sorted(missing_countries)}")

    return selected
