# SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Derives representative US transmission line types from the ACTIVSg82k network.

The script parses the PSS/E v33 RAW representation of the ACTIVSg82k synthetic
transmission network, converts branch electrical parameters into the units
required by PyPSA, and derives one representative transmission line type for
each supported nominal voltage level.

Inputs
------

- ``data/SyntheticUSA.RAW``:
    Native PSS/E v33 representation of the ACTIVSg82k synthetic transmission
    network.

Outputs
-------

- ``data/custom_line_types.csv``:
    Country-specific transmission line types, updated with the representative US transmission line types.

Description
-----------

Only in-service AC branches with positive length and positive RATEA are used.
Electrical parameters are converted from the PSS/E per-unit representation to
PyPSA line type units. The converted parameters are grouped by nominal voltage
and the median value is used to derive one representative line type per voltage
level.
"""

import csv
from pathlib import Path

import numpy as np
import pandas as pd

S_BASE_MVA = 100.0
MILE_TO_KM = 1.609344
FREQUENCY_HZ = 60.0
OMEGA = 2 * np.pi * FREQUENCY_HZ

COUNTRY = "US"
RAW_PATH = Path("data/SyntheticUSA.RAW")
OUTPUT_PATH = Path("data/custom_line_types.csv")

SUPPORTED_VOLTAGES = [69, 100, 115, 138, 161, 230, 345, 500, 765]


def read_raw_lines(raw_path: Path) -> list[str]:
    """
    Read a PSS/E RAW file as plain text lines.

    Parameters
    ----------
    raw_path : pathlib.Path
        Path to the PSS/E RAW file.

    Returns
    -------
    list[str]
        RAW file content split into individual lines.
    """
    if not raw_path.exists():
        msg = f"RAW file not found: {raw_path}"
        raise FileNotFoundError(msg)

    return raw_path.read_text(errors="ignore").splitlines()


def find_section_end(lines: list[str], marker: str) -> int:
    """
    Locate the end line of a RAW data section.

    Parameters
    ----------
    lines : list[str]
        RAW file lines.
    marker : str
        Section end marker, for example ``END OF BUS DATA``.

    Returns
    -------
    int
        Zero-based line index containing the requested marker.
    """
    try:
        return next(i for i, line in enumerate(lines) if marker in line)
    except StopIteration as exc:
        msg = f"Could not find RAW section marker: {marker}"
        raise ValueError(msg) from exc


def parse_buses(lines: list[str]) -> dict[int, float]:
    """
    Parse bus nominal voltages from the RAW bus section.

    Parameters
    ----------
    lines : list[str]
        RAW file lines.

    Returns
    -------
    dict[int, float]
        Mapping from bus ID to nominal voltage in kV.
    """
    bus_end = find_section_end(lines, "END OF BUS DATA")

    buses = {}
    for line in lines[1:bus_end]:
        raw = line.split("/")[0].strip()

        if not raw or not raw[0].isdigit():
            continue

        row = next(csv.reader([raw], skipinitialspace=True))
        bus_id = int(row[0])
        base_kv = float(row[2])
        buses[bus_id] = base_kv

    if not buses:
        msg = "No buses were parsed from the RAW file."
        raise ValueError(msg)

    return buses


def parse_branches(lines: list[str], buses: dict[int, float]) -> pd.DataFrame:
    """
    Parse AC branch records from the RAW branch section.

    Parameters
    ----------
    lines : list[str]
        RAW file lines.
    buses : dict[int, float]
        Mapping from bus ID to nominal voltage in kV.

    Returns
    -------
    pandas.DataFrame
        Branch table containing bus IDs, voltage levels, per-unit electrical
        parameters, thermal ratings, status, and length.
    """
    branch_start = find_section_end(lines, "END OF GENERATOR DATA") + 1
    branch_end = find_section_end(lines, "END OF BRANCH DATA")

    records = []

    for line in lines[branch_start:branch_end]:
        raw = line.split("/")[0].strip()

        if not raw:
            continue

        row = next(csv.reader([raw], skipinitialspace=True))

        if len(row) < 16:
            continue

        bus0 = int(str(row[0]).strip().strip("'\""))
        bus1 = int(str(row[1]).strip().strip("'\""))

        records.append(
            {
                "bus0": bus0,
                "bus1": bus1,
                "v_nom": buses[bus0],
                "v_nom_1": buses[bus1],
                "r_pu": float(row[3]),
                "x_pu": float(row[4]),
                "b_pu": float(row[5]),
                "rate_a": float(row[6]),
                "status": int(row[13]),
                "length_mile": float(row[15]),
            }
        )

    if not records:
        msg = "No AC branches were parsed from the RAW file."
        raise ValueError(msg)

    return pd.DataFrame.from_records(records)


def validate_branches(branches: pd.DataFrame) -> None:
    """
    Validate basic assumptions required by the derivation methodology.

    Parameters
    ----------
    branches : pandas.DataFrame
        Parsed branch table.

    Raises
    ------
    ValueError
        If any AC branch connects buses with different nominal voltage.
    """
    different_voltage = (branches["v_nom"] != branches["v_nom_1"]).sum()

    if different_voltage:
        msg = (
            f"Found {different_voltage} AC branches connecting different "
            "nominal voltage levels. Transformer records should not be part "
            "of the AC branch section."
        )
        raise ValueError(msg)


def convert_to_pypsa_units(branches: pd.DataFrame) -> pd.DataFrame:
    """
    Convert PSS/E branch parameters to PyPSA line type units.

    Parameters
    ----------
    branches : pandas.DataFrame
        Parsed branch table.

    Returns
    -------
    pandas.DataFrame
        Branch table extended with PyPSA-compatible electrical parameters:

        - ``r_per_length`` in Ohm/km,
        - ``x_per_length`` in Ohm/km,
        - ``c_per_length`` in nF/km,
        - ``i_nom`` in kA.
    """
    branches = branches.copy()

    branches["length_km"] = branches["length_mile"] * MILE_TO_KM

    z_base = branches["v_nom"] ** 2 / S_BASE_MVA
    y_base = S_BASE_MVA / branches["v_nom"] ** 2

    branches["r_per_length"] = branches["r_pu"] * z_base / branches["length_km"]
    branches["x_per_length"] = branches["x_pu"] * z_base / branches["length_km"]

    b_siemens = branches["b_pu"] * y_base
    c_farads = b_siemens / OMEGA

    branches["c_per_length"] = c_farads / branches["length_km"] * 1e9
    branches["i_nom"] = branches["rate_a"] / (np.sqrt(3) * branches["v_nom"])

    return branches


def filter_valid_branches(branches: pd.DataFrame) -> pd.DataFrame:
    """
    Filter branch records used for representative line type derivation.

    Parameters
    ----------
    branches : pandas.DataFrame
        Branch table after unit conversion.

    Returns
    -------
    pandas.DataFrame
        Valid branch records with active status, positive length, positive
        RATEA, equal terminal voltages, and supported nominal voltage.
    """
    valid = branches[
        (branches["status"] == 1)
        & (branches["length_km"] > 0)
        & (branches["rate_a"] > 0)
        & (branches["v_nom"] == branches["v_nom_1"])
        & (branches["v_nom"].isin(SUPPORTED_VOLTAGES))
    ].copy()

    if valid.empty:
        msg = "No valid branch records remain after filtering."
        raise ValueError(msg)

    return valid


def build_line_types(branches: pd.DataFrame) -> pd.DataFrame:
    """
    Derive one representative line type per supported voltage level.

    Parameters
    ----------
    branches : pandas.DataFrame
        Branch table after conversion to PyPSA units.

    Returns
    -------
    pandas.DataFrame
        Representative line types indexed by stable type names.
    """
    valid = filter_valid_branches(branches)

    line_types = (
        valid.groupby("v_nom")
        .agg(
            r_per_length=("r_per_length", "median"),
            x_per_length=("x_per_length", "median"),
            c_per_length=("c_per_length", "median"),
            i_nom=("i_nom", "median"),
            samples=("v_nom", "size"),
        )
        .reset_index()
    )

    missing = sorted(set(SUPPORTED_VOLTAGES) - set(line_types["v_nom"].astype(int)))

    if missing:
        msg = f"Missing representative line types for voltages: {missing}"
        raise ValueError(msg)

    line_types["type"] = line_types["v_nom"].astype(int).map(lambda v: f"US_{v}kV")
    line_types["f_nom"] = FREQUENCY_HZ

    line_types = line_types[
        [
            "type",
            "v_nom",
            "r_per_length",
            "x_per_length",
            "c_per_length",
            "i_nom",
            "f_nom",
            "samples",
        ]
    ]

    return line_types.sort_values("v_nom")


def write_line_types(line_types: pd.DataFrame, output_path: Path) -> None:
    """
    Update the US definitions in the custom line-type registry.

    Existing definitions for other countries are preserved, while previous US
    definitions are replaced with the newly derived values.

    Parameters
    ----------
    line_types : pandas.DataFrame
        Representative US line type table.
    output_path : pathlib.Path
        Destination custom line-type registry.
    """
    line_types = line_types.copy()
    line_types.insert(0, "country", COUNTRY)

    if output_path.exists():
        existing = pd.read_csv(output_path)

        required_registry_columns = {"country", "type"}
        missing_columns = required_registry_columns - set(existing.columns)

        if missing_columns:
            raise ValueError(
                f"Custom line type registry {output_path} is missing columns: "
                f"{sorted(missing_columns)}"
            )

        existing = existing.loc[existing["country"] != COUNTRY]
        line_types = pd.concat([existing, line_types], ignore_index=True)

    duplicated_types = line_types.loc[line_types["type"].duplicated(), "type"].unique()

    if len(duplicated_types):
        raise ValueError(
            f"Duplicate line types in custom registry: " f"{duplicated_types.tolist()}"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    line_types = line_types.sort_values(
        ["country", "v_nom"],
        kind="stable",
    )

    line_types.to_csv(output_path, index=False)


def main() -> None:
    """
    Run the full ACTIVSg82k line type derivation workflow.
    """
    lines = read_raw_lines(RAW_PATH)
    buses = parse_buses(lines)
    branches = parse_branches(lines, buses)

    validate_branches(branches)

    branches = convert_to_pypsa_units(branches)
    line_types = build_line_types(branches)

    write_line_types(line_types, OUTPUT_PATH)

    print(line_types.to_string(index=False))
    print(f"\nWritten: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
