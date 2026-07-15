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


def make_grouped_bar_html(carriers_data, unit):
    """
    Generate clustered horizontal bar chart in HTML.
    carriers_data: {carrier -> {"model": float, "sources": {source -> float}}}
    unit: 'MW' or 'TWh'
    """
    if not carriers_data:
        return ""

    # Find max value to scale the bars
    max_val = 0.0
    for data in carriers_data.values():
        max_val = max(max_val, data["model"])
        for val in data["sources"].values():
            max_val = max(max_val, val)

    if max_val <= 0:
        return ""

    # Sort carriers by model value descending
    sorted_carriers = sorted(
        carriers_data.keys(), key=lambda c: carriers_data[c]["model"], reverse=True
    )

    # Color mapping for sources
    source_colors = {
        "model": "#2563eb",
        "ember": "#10b981",
        "irena": "#f59e0b",
        "ourworldindata": "#8b5cf6",
        "owid": "#8b5cf6",
    }

    rows_html = []
    for c in sorted_carriers:
        data = carriers_data[c]

        # Skip displaying this carrier if both model and reference values are all zero
        model_val = data["model"]
        ref_vals = data["sources"].values()
        if model_val <= 0 and all(v <= 0 for v in ref_vals):
            continue

        carrier_name = c.upper()

        # Build list of bars to render for this carrier
        bars = [("Model", data["model"], source_colors["model"])]
        for src, val in sorted(data["sources"].items()):
            color = source_colors.get(src.lower(), "#64748b")
            label = "OWID" if src.lower() == "ourworldindata" else src.upper()
            bars.append((label, val, color))

        # Render bars for this carrier
        bars_html = []
        for src_label, val, color in bars:
            width_pct = (val / max_val) * 100 if max_val > 0 else 0
            val_str = f"{val:.1f}" if val > 0 else "0"

            bars_html.append(
                f"""
            <div style="display: flex; align-items: center; margin-bottom: 3px;">
                <span style="font-size: 10px; font-weight: 500; color: #64748b; width: 60px; text-transform: uppercase;">{src_label}</span>
                <div style="flex-grow: 1; background-color: #f1f5f9; height: 10px; border-radius: 5px; overflow: hidden; border: 1px solid #e2e8f0; max-width: 350px;">
                    <div style="width: {width_pct:.1f}%; background-color: {color}; height: 100%; border-radius: 4px;" title="{src_label}: {val_str} {unit}"></div>
                </div>
                <span style="font-size: 10px; font-weight: 600; color: #334155; margin-left: 8px; width: 80px;">{val_str} {unit}</span>
            </div>
            """
            )

        rows_html.append(
            f"""
        <div style="display: flex; margin-bottom: 16px; border-bottom: 1px solid #f1f5f9; padding-bottom: 12px; align-items: flex-start;">
            <div style="width: 110px; font-weight: 600; color: #475569; font-size: 12px; padding-top: 2px; text-align: right; margin-right: 16px; word-break: break-all;">
                {carrier_name}
            </div>
            <div style="flex-grow: 1;">
                {"".join(bars_html)}
            </div>
        </div>
        """
        )

    # Legend HTML
    legend_items = []
    for label, color in [
        ("Model", "#2563eb"),
        ("EMBER", "#10b981"),
        ("IRENA", "#f59e0b"),
        ("OWID", "#8b5cf6"),
    ]:
        legend_items.append(
            f"""
        <span style="display: inline-flex; align-items: center; margin-right: 16px; font-size: 11px; color: #475569; font-weight: 500;">
            <i style="display: inline-block; width: 12px; height: 12px; border-radius: 3px; background-color: {color}; margin-right: 6px;"></i>
            {label}
        </span>
        """
        )

    html = f"""
<div style="margin: 20px 0; padding: 16px; border: 1px solid #e2e8f0; border-radius: 8px; background-color: #f8fafc; max-width: 650px; box-shadow: 0 1px 3px rgba(0,0,0,0.02);">
    <div style="display: flex; flex-wrap: wrap; margin-bottom: 16px; padding-bottom: 10px; border-bottom: 2px solid #e2e8f0;">
        {"".join(legend_items)}
    </div>
    {"".join(rows_html)}
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
            demand_carriers = {}
            for r in grp:
                if r["pillar"] == "demand":
                    py_val = parse_float(r["pypsa_value"]) or 0.0
                    ref_val = parse_float(r["reference_value"]) or 0.0
                    source = r["reference_source"]
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
                    if "total demand" not in demand_carriers:
                        demand_carriers["total demand"] = {
                            "model": py_val,
                            "sources": {},
                        }
                    demand_carriers["total demand"]["sources"][source] = ref_val

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

                chart_html = make_grouped_bar_html(demand_carriers, "TWh")
                if chart_html:
                    f.write("**Total Demand Comparison:**\n\n")
                    f.write(chart_html)
                    f.write("\n")

            # Capacity section
            capacity_rows = []
            cap_carriers = {}
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
                        if carrier not in cap_carriers:
                            cap_carriers[carrier] = {"model": py_val, "sources": {}}
                        cap_carriers[carrier]["sources"][source] = ref_val

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

                chart_html = make_grouped_bar_html(cap_carriers, "MW")
                if chart_html:
                    f.write("**Installed Capacity Mix Comparison:**\n\n")
                    f.write(chart_html)
                    f.write("\n")

            # Generation section
            generation_rows = []
            gen_carriers = {}
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
                        if carrier not in gen_carriers:
                            gen_carriers[carrier] = {"model": py_val, "sources": {}}
                        gen_carriers[carrier]["sources"][source] = ref_val

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

                chart_html = make_grouped_bar_html(gen_carriers, "TWh")
                if chart_html:
                    f.write("**Electricity Generation Mix Comparison:**\n\n")
                    f.write(chart_html)
                    f.write("\n")

            f.write("---\n\n")

    logger.info("Validation Dashboard Hook: Successfully completed page compilation.")
