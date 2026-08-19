# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Render configs/datasources_url_map.toml into a human-readable markdown page.

Writes doc/user-guide/data_api.md, with one second-level section per
dataset (titled with its "long_name"), starting with the local "output"
path and followed by the "description".

Usage
-----
python scripts/non_workflow/generate_data_api_doc.py
"""

import tomllib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_TOML = REPO_ROOT / "configs" / "datasources_url_map.toml"
OUTPUT_MD = REPO_ROOT / "doc" / "user-guide" / "data_api.md"

# Tutorial-scoped datasets (entries with tutorial = true in the toml) are
# always skipped -- they're bundle-specific copies of the entries below.

# Order in which datasets are rendered. Edit this list by hand to reorder
# the doc; any dataset name present in the toml but missing here is
# appended at the end (in its original toml order) rather than dropped.
DATASET_ORDER = [
    "gadm",
    "gadm_v36",
    "worldpop_maxar",
    "worldpop",
    "worldpop_api",
    "demandcast_forecasts",
    "gegis_demand_projections",
    "eez_marineregions",
    "gebco_bathymetry",
    "copernicus_landcover",
    "wdpa_protectedplanet",
    "natura_raster",
    "hydrobasins",
    "edgar",
    "irena_statistics",
    "global_buildings_microsoft",
    "global_buildings_microsoft_quadrants",
    "un_energy_balances_unsd",
    "pop_total_un",
    "airports",
    "air_runways",
    "steel_gem",
    "pipelines_gem",
    "gas_network_iggielgn",
    "refineries",
    "osm_nominatim_geocoding",
    "n_vehicles_who",
    "vehicles_per_capita_wiki",
    "transport_emission_worldbank",
    "sea_ports_nga",
]

MD_HEADER = """<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

# Description of datasets used by workflow

"""


def render_entry(entry):
    lines = [f"## {entry['long_name']}", ""]

    output = entry.get("output")
    if output:
        lines.append(f"**Output:** `{output}`")
        lines.append("")

    description = entry.get("description")
    if description:
        lines.append(description)
        lines.append("")

    return "\n".join(lines)


def order_entries(entries):
    by_name = {entry["name"]: entry for entry in entries}
    ordered = [by_name.pop(name) for name in DATASET_ORDER if name in by_name]
    # anything not listed in DATASET_ORDER is appended, in its original order,
    # instead of being silently dropped
    ordered.extend(by_name.values())
    return ordered


def main():
    with open(SOURCE_TOML, "rb") as f:
        data = tomllib.load(f)

    entries = [entry for entry in data["source"] if not entry.get("tutorial")]
    entries = order_entries(entries)
    sections = [render_entry(entry) for entry in entries]

    OUTPUT_MD.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_MD, "w") as f:
        f.write(MD_HEADER)
        f.write("\n".join(sections))

    print(f"Wrote {len(sections)} dataset sections to {OUTPUT_MD}")


if __name__ == "__main__":
    main()
