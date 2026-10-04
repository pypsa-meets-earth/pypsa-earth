# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Render configs/datasources_url_map.yaml into a description of the datasets
in the documentation.

Usage
-----
python doc/assets/scripts/generate_data_api_doc.py
"""

import logging
from pathlib import Path

import yaml

logger = logging.getLogger(__name__)

SOURCE_YAML = Path("configs/datasources_url_map.yaml")
OUTPUT_MD = Path("doc/user-guide/data_api.md")

# Tutorial-scoped datasets (entries with tutorial: true in the yaml) are
# always skipped -- they're bundle-specific copies of the entries below.

# Order in which datasets are rendered. Edit this list by hand to reorder
# the doc; any dataset name present in the yaml but missing here is
# appended at the end (in its original yaml order) rather than dropped.
DATASET_ORDER = [
    "osm_geofabrik",
    "era5",
    "sarah3",
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

# Description of datasets used by PyPSA-Earth workflow

"""


TABLE_HEADER = "| Dataset | Output | Description |\n|---|---|---|\n"


def render_cell(value):
    return str(value).replace("|", "\\|").replace("\n", " ")


def render_row(entry):
    output = entry.get("output")
    output = f"`{render_cell(output)}`" if output else ""
    description = render_cell(entry.get("description", ""))
    return f"| {render_cell(entry['long_name'])} | {output} | {description} |"


def order_entries(entries, order=DATASET_ORDER):
    by_name = {entry["name"]: entry for entry in entries}
    ordered = [by_name.pop(name) for name in order if name in by_name]
    ordered.extend(by_name.values())
    return ordered


def main():
    with open(SOURCE_YAML) as f:
        data = yaml.safe_load(f)

    entries = [entry for entry in data["source"] if not entry.get("tutorial")]
    entries = order_entries(entries)
    rows = [render_row(entry) for entry in entries]

    OUTPUT_MD.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_MD, "w") as f:
        f.write(MD_HEADER)
        f.write(TABLE_HEADER)
        f.write("\n".join(rows) + "\n")

    logger.info(f"Wrote {len(rows)} dataset rows to {OUTPUT_MD}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
