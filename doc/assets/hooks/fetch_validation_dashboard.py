# SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
MkDocs hook to automatically download validation results from the
pypsa-earth-status repository and generate the validation dashboard pages.
"""

import csv
import logging
import os
import urllib.request
from pathlib import Path

logger = logging.getLogger("mkdocs")

# TODO: Update to production repository URL after the pypsa-earth-status PR is merged.
CSV_URL = (
    "https://raw.githubusercontent.com/GbotemiB/pypsa-earth-status/"
    "health-status/results/health_status.csv"
)

# World boundaries GeoJSON for Leaflet
GEOJSON_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/"
    "geojson/ne_110m_admin_0_countries.geojson"
)


def make_markdown_table(headers, alignments, rows):
    """Generate a markdown table with custom alignment."""
    separator = []
    for align in alignments:
        if align == "right":
            separator.append("---:")
        elif align == "center":
            separator.append(":---:")
        else:
            separator.append(":---")
    lines = []
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join(separator) + " |")
    for r in rows:
        lines.append("| " + " | ".join(str(r.get(h, "")) for h in headers) + " |")
    return "\n".join(lines) + "\n"


def parse_float(val):
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


CARRIER_COLORS = {
    "coal": "#374151",
    "ccgt": "#f97316",
    "gas": "#f97316",
    "oil": "#dc2626",
    "nuclear": "#c084fc",
    "pv": "#facc15",
    "solar": "#facc15",
    "wind": "#3b82f6",
    "hydro": "#06b6d4",
    "biomass": "#22c55e",
    "geothermal": "#b45309",
    "other": "#94a3b8",
}


def make_stacked_bar_html(mix_dict, unit):
    """
    Generate horizontal stacked bar charts in HTML for markdown page comparison.
    mix_dict: {carrier -> (pypsa_val, ref_val)}
    unit: 'MW' or 'TWh'
    """
    total_py = sum(val[0] for val in mix_dict.values())
    total_ref = sum(val[1] for val in mix_dict.values())

    if total_py <= 0 and total_ref <= 0:
        return ""

    # Sort carriers by reference value descending
    sorted_carriers = sorted(
        mix_dict.keys(), key=lambda c: mix_dict[c][1], reverse=True
    )

    # Render segments
    segments_py = []
    segments_ref = []
    legend_items = []

    for c in sorted_carriers:
        py_val, ref_val = mix_dict[c]
        color = CARRIER_COLORS.get(c.lower(), "#94a3b8")

        # PyPSA Segment
        if total_py > 0 and py_val > 0:
            pct_py = (py_val / total_py) * 100
            segments_py.append(
                f'<div style="width: {pct_py:.1f}%; background-color: {color}; height: 100%;" '
                f'title="{c.upper()}: {py_val:.1f} {unit} ({pct_py:.1f}%)"></div>'
            )

        # Ref Segment
        if total_ref > 0 and ref_val > 0:
            pct_ref = (ref_val / total_ref) * 100
            segments_ref.append(
                f'<div style="width: {pct_ref:.1f}%; background-color: {color}; height: 100%;" '
                f'title="{c.upper()}: {ref_val:.1f} {unit} ({pct_ref:.1f}%)"></div>'
            )

        # Legend
        if py_val > 0 or ref_val > 0:
            legend_items.append(
                f'<span style="display: inline-flex; align-items: center; margin-right: 12px; font-size: 11px; color: #475569;">'
                f'<i style="display: inline-block; width: 10px; height: 10px; border-radius: 2px; background-color: {color}; margin-right: 4px;"></i>'
                f"{c.upper()}"
                f"</span>"
            )

    py_bars = "".join(segments_py)
    ref_bars = "".join(segments_ref)
    legend_html = "".join(legend_items)

    html = f"""
<div style="margin: 16px 0; padding: 12px; border: 1px solid #e2e8f0; border-radius: 6px; background-color: #f8fafc; max-width: 600px;">
    <div style="display: flex; align-items: center; margin-bottom: 6px;">
        <span style="font-size: 11px; font-weight: 600; color: #64748b; width: 80px;">Model:</span>
        <div style="flex-grow: 1; display: flex; height: 16px; border-radius: 3px; overflow: hidden; background-color: #e2e8f0; border: 1px solid #cbd5e1;">
            {py_bars or '<div style="width: 100%; background-color: #e2e8f0; height: 100%; font-size: 9px; text-align: center; color: #94a3b8; line-height: 14px;">No Generation</div>'}
        </div>
        <span style="font-size: 11px; font-weight: 600; color: #334155; margin-left: 8px; width: 70px; text-align: right;">{total_py:.1f} {unit}</span>
    </div>
    <div style="display: flex; align-items: center; margin-bottom: 8px;">
        <span style="font-size: 11px; font-weight: 600; color: #64748b; width: 80px;">Reference:</span>
        <div style="flex-grow: 1; display: flex; height: 16px; border-radius: 3px; overflow: hidden; background-color: #e2e8f0; border: 1px solid #cbd5e1;">
            {ref_bars or '<div style="width: 100%; background-color: #e2e8f0; height: 100%; font-size: 9px; text-align: center; color: #94a3b8; line-height: 14px;">No Generation</div>'}
        </div>
        <span style="font-size: 11px; font-weight: 600; color: #334155; margin-left: 8px; width: 70px; text-align: right;">{total_ref:.1f} {unit}</span>
    </div>
    <div style="display: flex; flex-wrap: wrap; margin-top: 8px; border-top: 1px solid #e2e8f0; padding-top: 6px;">
        {legend_html}
    </div>
</div>
"""
    return html


MAP_TEMPLATE = ""


def on_pre_build(config, **kwargs):
    """MkDocs hook: fetch validation data and generate dashboard pages."""
    logger.info("Validation Dashboard Hook: Starting build setup...")

    # Define paths
    docs_dir = Path(config["docs_dir"])
    val_dir = docs_dir / "validation"
    val_dir.mkdir(parents=True, exist_ok=True)

    csv_path = val_dir / "health_status.csv"
    overview_path = val_dir / "overview.md"
    statistics_path = val_dir / "statistics.md"

    # 1. Fetch CSV (only if not already present locally)
    if not csv_path.exists():
        try:
            logger.info(
                f"Validation Dashboard Hook: Fetching validation data from {CSV_URL}"
            )
            with urllib.request.urlopen(CSV_URL, timeout=10) as response:
                csv_data = response.read().decode("utf-8")
            with open(csv_path, "w", encoding="utf-8") as f:
                f.write(csv_data)
            logger.info(
                "Validation Dashboard Hook: Successfully saved health_status.csv"
            )
        except Exception as e:
            logger.warning(
                f"Validation Dashboard Hook: Failed to download CSV ({e}). "
                "Falling back to existing local copy if available."
            )
            if not csv_path.exists():
                logger.error(
                    "Validation Dashboard Hook: No fallback CSV found! Skipping page generation."
                )
                return
    else:
        logger.info(
            "Validation Dashboard Hook: Using existing local copy of health_status.csv"
        )

    # 2. Parse CSV
    records = []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            records.append(row)

    if not records:
        logger.warning(
            "Validation Dashboard Hook: CSV file is empty! Skipping page generation."
        )
        return

    # Group by scenario key / country details
    groups = {}
    for r in records:
        key = (
            r["scenario_key"],
            r["country_code"],
            r["country_name"],
            r["pypsa_earth_version"],
            r["year"],
        )
        if key not in groups:
            groups[key] = []
        groups[key].append(r)

    # Sort groups by country name
    sorted_keys = sorted(groups.keys(), key=lambda k: k[2])

    # 3. Generate overview.md
    logger.info("Validation Dashboard Hook: Generating overview.md...")
    overview_rows = []
    for key in sorted_keys:
        scenario, country_code, country_name, version, year = key
        grp = groups[key]

        def format_grades(pillar, metric):
            grades = []
            for r in grp:
                if r["pillar"] == pillar and r["metric"] == metric:
                    grade = r["grade"]
                    src = r["reference_source"]
                    if grade and grade.strip():
                        grades.append(f"{grade} ({src})")
            return ", ".join(grades) if grades else "-"

        anchor_slug = f"{country_name.lower().replace(' ', '-')}-{country_code.lower()}-scenario-{scenario.lower()}"
        overview_rows.append(
            {
                "Scenario": scenario,
                "Country": f"[{country_name} ({country_code})](statistics.md#{anchor_slug})",
                "Demand Grade": format_grades("demand", "total_demand"),
                "Capacity Grade": format_grades("installed_capacity", "total_capacity"),
                "Generation Grade": format_grades("generation", "total_generation"),
                "PyPSA-Earth Version": version,
                "Year": year,
            }
        )

    headers = [
        "Scenario",
        "Country",
        "Demand Grade",
        "Capacity Grade",
        "Generation Grade",
        "PyPSA-Earth Version",
        "Year",
    ]
    alignments = ["left", "left", "left", "left", "left", "left", "right"]

    with open(overview_path, "w", encoding="utf-8") as f:
        f.write("# Validation Overview\n\n")
        f.write(
            "A summary of validation grades across all analyzed countries and scenarios. "
            "Grades assess the percentage deviation between model outputs and historical reference statistics (Ember, IRENA, OWID).\n\n"
        )
        f.write(make_markdown_table(headers, alignments, overview_rows))

    # 4. Generate statistics.md
    logger.info("Validation Dashboard Hook: Generating statistics.md...")
    with open(statistics_path, "w", encoding="utf-8") as f:
        f.write("# Detailed Validation Statistics\n\n")
        f.write(
            "Detailed comparisons of active validation scenarios. Metrics are grouped by country and pillar.\n\n"
        )

        for key in sorted_keys:
            scenario, country_code, country_name, version, year = key
            grp = groups[key]

            f.write(f"## {country_name} ({country_code}) — Scenario `{scenario}`\n\n")
            f.write(f"* **Model Year:** {year}\n")
            f.write(f"* **PyPSA-Earth Version:** `{version}`\n\n")

            # Demand section
            demand_rows = []
            for r in grp:
                if r["pillar"] == "demand":
                    py_val = parse_float(r["pypsa_value"])
                    ref_val = parse_float(r["reference_value"])
                    dev = parse_float(r["deviation_pct"])
                    dev_str = f"{dev:+.2f}%" if dev is not None else "-"

                    demand_rows.append(
                        {
                            "Source": r["reference_source"],
                            "Model Value (TWh)": (
                                f"{py_val:.2f}" if py_val is not None else "-"
                            ),
                            "Reference Value (TWh)": (
                                f"{ref_val:.2f}" if ref_val is not None else "-"
                            ),
                            "Deviation (%)": dev_str,
                            "Grade": r["grade"] if r["grade"] else "-",
                        }
                    )

            if demand_rows:
                f.write("### 1. Electricity Demand\n\n")
                f.write(
                    make_markdown_table(
                        [
                            "Source",
                            "Model Value (TWh)",
                            "Reference Value (TWh)",
                            "Deviation (%)",
                            "Grade",
                        ],
                        ["left", "right", "right", "left", "left"],
                        demand_rows,
                    )
                )
                f.write("\n")

            # Capacity section
            capacity_rows = []
            cap_mix = {}
            for r in grp:
                if r["pillar"] == "installed_capacity":
                    if r["metric"] == "total_capacity":
                        py_val = parse_float(r["pypsa_value"])
                        ref_val = parse_float(r["reference_value"])
                        dev = parse_float(r["deviation_pct"])
                        dev_str = f"{dev:+.2f}%" if dev is not None else "-"

                        capacity_rows.append(
                            {
                                "Source": r["reference_source"],
                                "Model Value (MW)": (
                                    f"{py_val:.2f}" if py_val is not None else "-"
                                ),
                                "Reference Value (MW)": (
                                    f"{ref_val:.2f}" if ref_val is not None else "-"
                                ),
                                "Deviation (%)": dev_str,
                                "Grade": r["grade"] if r["grade"] else "-",
                            }
                        )
                    elif r["metric"] == "capacity":
                        source = r["reference_source"]
                        carrier = r["carrier"]
                        py_val = parse_float(r["pypsa_value"]) or 0.0
                        ref_val = parse_float(r["reference_value"]) or 0.0
                        if source not in cap_mix:
                            cap_mix[source] = {}
                        cap_mix[source][carrier] = (py_val, ref_val)

            if capacity_rows:
                f.write("### 2. Installed Capacity\n\n")
                f.write("**Total Installed Capacity Comparison:**\n\n")
                f.write(
                    make_markdown_table(
                        [
                            "Source",
                            "Model Value (MW)",
                            "Reference Value (MW)",
                            "Deviation (%)",
                            "Grade",
                        ],
                        ["left", "right", "right", "left", "left"],
                        capacity_rows,
                    )
                )
                f.write("\n")

                for source, mix in cap_mix.items():
                    chart_html = make_stacked_bar_html(mix, "MW")
                    if chart_html:
                        f.write(f"**Installed Capacity Mix ({source.upper()}):**\n\n")
                        f.write(chart_html)
                        f.write("\n")

            # Generation section
            generation_rows = []
            gen_mix = {}
            for r in grp:
                if r["pillar"] == "generation":
                    if r["metric"] == "total_generation":
                        py_val = parse_float(r["pypsa_value"])
                        ref_val = parse_float(r["reference_value"])
                        dev = parse_float(r["deviation_pct"])
                        dev_str = f"{dev:+.2f}%" if dev is not None else "-"

                        generation_rows.append(
                            {
                                "Source": r["reference_source"],
                                "Model Value (TWh)": (
                                    f"{py_val:.2f}" if py_val is not None else "-"
                                ),
                                "Reference Value (TWh)": (
                                    f"{ref_val:.2f}" if ref_val is not None else "-"
                                ),
                                "Deviation (%)": dev_str,
                                "Grade": r["grade"] if r["grade"] else "-",
                            }
                        )
                    elif r["metric"] == "generation":
                        source = r["reference_source"]
                        carrier = r["carrier"]
                        py_val = parse_float(r["pypsa_value"]) or 0.0
                        ref_val = parse_float(r["reference_value"]) or 0.0
                        if source not in gen_mix:
                            gen_mix[source] = {}
                        gen_mix[source][carrier] = (py_val, ref_val)

            if generation_rows:
                f.write("### 3. Electricity Generation\n\n")
                f.write("**Total Generation Comparison:**\n\n")
                f.write(
                    make_markdown_table(
                        [
                            "Source",
                            "Model Value (TWh)",
                            "Reference Value (TWh)",
                            "Deviation (%)",
                            "Grade",
                        ],
                        ["left", "right", "right", "left", "left"],
                        generation_rows,
                    )
                )
                f.write("\n")

                for source, mix in gen_mix.items():
                    chart_html = make_stacked_bar_html(mix, "TWh")
                    if chart_html:
                        f.write(
                            f"**Electricity Generation Mix ({source.upper()}):**\n\n"
                        )
                        f.write(chart_html)
                        f.write("\n")

            f.write("---\n\n")

    logger.info("Validation Dashboard Hook: Successfully completed page compilation.")
