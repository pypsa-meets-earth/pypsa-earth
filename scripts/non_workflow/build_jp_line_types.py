# SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Build representative Japanese transmission line types.

Japan uses transmission voltage classes that differ from the default
European line-type mapping used by PyPSA-Earth. This utility creates one
representative line type for each major Japanese voltage class and updates
only the Japanese entries in the shared custom line-type registry.

The voltage classes are based on public Japanese transmission-system
documentation from OCCTO and regional transmission utilities.

The electrical parameters are representative engineering defaults rather
than line-specific measurements. The reference values are based on the
open All-Japan-Grid transmission-line table and are checked for consistency
with publicly available Japanese engineering references. In particular,
the 500 kV values are consistent with the IEEJ standard bulk-power-system
model using a four-conductor TACSR 810 mm2 transmission line.

Japan operates both 50 Hz and 60 Hz AC systems. PyPSA-Earth currently
selects country-specific line types by country and voltage only and does
not preserve line frequency through the complete OSM-to-PyPSA pipeline.
The initial Japanese line types therefore use 50 Hz as the reference
frequency when converting shunt susceptance to capacitance.

This is a non-workflow utility. Running it replaces only JP entries in
``data/custom_line_types.csv`` and preserves definitions for other
countries unchanged.

Sources
-------
IEEJ bulk power system standard model:
https://www.iee.or.jp/pes/model/english/kikan/Models/index.html

OCCTO grid information:
https://www.occto.or.jp/en/information_disclosure/

All-Japan-Grid:
https://github.com/lutelute/All-Japan-Grid
"""

import csv
from pathlib import Path

import numpy as np

COUNTRY = "JP"
FREQUENCY_HZ = 50.0
OUTPUT_PATH = Path("data/custom_line_types.csv")

COLUMNS = [
    "country",
    "type",
    "v_nom",
    "r_per_length",
    "x_per_length",
    "c_per_length",
    "i_nom",
    "f_nom",
    "samples",
]

# v_nom [kV], r/x [Ohm/km], b [S/km], i_nom [kA]
REFERENCE_TYPES = [
    (66, 0.120, 0.400, 3.20e-6, 0.6),
    (77, 0.100, 0.395, 3.30e-6, 0.7),
    (110, 0.055, 0.385, 3.45e-6, 0.9),
    (132, 0.045, 0.370, 3.55e-6, 1.2),
    (154, 0.050, 0.380, 3.50e-6, 1.0),
    (187, 0.038, 0.350, 3.65e-6, 1.5),
    (220, 0.032, 0.335, 3.75e-6, 1.8),
    (275, 0.028, 0.325, 3.85e-6, 2.0),
    (500, 0.012, 0.290, 4.10e-6, 4.0),
]


def build_line_types() -> list[list[str]]:
    """Build representative Japanese line types in PyPSA registry units."""
    omega = 2 * np.pi * FREQUENCY_HZ
    rows = []

    for v_nom, resistance, reactance, susceptance, current in REFERENCE_TYPES:
        if min(resistance, reactance, susceptance, current) <= 0:
            raise ValueError(f"Line-type parameters for {v_nom} kV must be positive.")

        capacitance = susceptance / omega * 1e9

        rows.append(
            [
                COUNTRY,
                f"JP_{v_nom}kV",
                str(v_nom),
                str(resistance),
                str(reactance),
                str(capacitance),
                str(current),
                str(int(FREQUENCY_HZ)),
                "",
            ]
        )

    return rows


def read_registry(path: Path) -> tuple[str, list[str]]:
    """Read the registry while preserving non-Japanese rows verbatim."""
    if not path.exists():
        raise FileNotFoundError(f"Custom line-type registry not found: {path}")

    lines = path.read_text().splitlines()

    if not lines:
        raise ValueError(f"Custom line-type registry is empty: {path}")

    header = lines[0]

    if next(csv.reader([header])) != COLUMNS:
        raise ValueError(f"Unexpected columns in {path}: {next(csv.reader([header]))}")

    other_rows = []

    for line in lines[1:]:
        if not line.strip():
            continue

        row = next(csv.reader([line]))

        if len(row) != len(COLUMNS):
            raise ValueError(f"Malformed registry row: {line}")

        if row[0] != COUNTRY:
            other_rows.append(line)

    return header, other_rows


def write_line_types(path: Path, rows: list[list[str]]) -> None:
    """Replace JP entries while leaving other country definitions unchanged."""
    header, other_rows = read_registry(path)

    jp_lines = [",".join(row) for row in rows]

    duplicate_types = [
        row[1]
        for index, row in enumerate(rows)
        if row[1] in {previous[1] for previous in rows[:index]}
    ]

    if duplicate_types:
        raise ValueError(f"Duplicate Japanese line types: {duplicate_types}")

    output = [header, *jp_lines, *other_rows]
    path.write_text("\n".join(output) + "\n")


def main() -> None:
    """Generate and write Japanese representative line types."""
    rows = build_line_types()
    write_line_types(OUTPUT_PATH, rows)

    print("\n".join(",".join(row) for row in rows))
    print(f"\nWritten: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
