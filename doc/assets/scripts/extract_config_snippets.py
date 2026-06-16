#!/usr/bin/env python3
"""
Extract YAML sections from config.default.yaml and save them as individual snippet files.
These snippets are embedded in the documentation (e.g., configuration.md) via pymdownx.snippets.

When to edit this script:
    - When a new top-level section is added to config.default.yaml
    - When a config section is renamed or removed
    - When a new subsection needs its own snippet file
    (Note: This script should be updated when modifying the config structure, as guided by the PR template.)

How it runs:
    - Automatically via the 'extract-config-snippets' pre-commit hook
    - Manually: python doc/assets/scripts/extract_config_snippets.py
"""

import re
from pathlib import Path


def extract_yaml_section(yaml_content, section_key, subsection=None):
    """
    Extract a section from YAML content preserving formatting and comments.
    Returns the text of that section.
    """
    lines = yaml_content.split("\n")
    section_lines = []
    in_section = False
    section_indent = None

    target_key = subsection if subsection else section_key

    for line in lines:
        # Check if we're starting the target section
        if re.match(rf"^{re.escape(target_key)}:", line.strip()):
            in_section = True
            section_indent = len(line) - len(line.lstrip())
            section_lines.append(line)
            continue

        if in_section:
            # Check if we've hit a new top-level or same-level key
            if line.strip() and not line.strip().startswith("#"):
                current_indent = len(line) - len(line.lstrip())
                # If we're back to the same or lower indent level, we're done
                if current_indent <= section_indent:
                    break

            section_lines.append(line)

    return "\n".join(section_lines)


def main():
    config_file = Path("config.default.yaml")
    snippets_dir = Path("doc/configtables/snippets")
    snippets_dir.mkdir(exist_ok=True)

    # Read the config file
    with open(config_file, "r") as f:
        config_content = f.read()

    # Mapping of snippet files to config sections
    sections = {
        "meta": {
            "start": "version:",
            "end": "# =================== STUDY SETUP ===================",
        },
        "study_setup": {
            "start": "countries:",
            "end": "run:",
        },
        "data_retrieval": {
            "start": "enable:",
            "end": "# =================== GEOGRAPHY & SHAPES ===================",
        },
        "run": ["run"],
        "scenario": ["scenario"],
        "snapshots": {
            "start": "snapshots:",
            "end": "# =================== DATA RETRIEVAL ===================",
        },
        "crs": {
            "start": "crs:",
            "end": "# ------------------- Regional shapes",
        },
        "augmented_line_connection": {
            "start": "augmented_line_connection:",
            "end": "# =================== NETWORK & RESOURCES",
        },
        "cluster_options": {
            "start": "cluster_options:",
            "end": "# ------------------- Augmented connectivity",
        },
        "build_shape_options": {
            "start": "build_shape_options:",
            "end": "# ------------------- Subregions",
        },
        "subregion": {
            "start": "subregion:",
            "end": "# ------------------- Land-cover exclusions",
        },
        "natura": {
            "start": "natura:",
            "end": "# ------------------- OpenStreetMap",
        },
        "osm": {
            "start": "osm:",
            "end": "# ------------------- Clustering",
        },
        "base_network": {
            "start": "base_network:",
            "end": "# ------------------- Demand",
        },
        "load_options": {
            "start": "load_options:",
            "end": "# ------------------- Electricity grid",
        },
        "electricity": {
            "start": "electricity:",
            "end": "lines:",
        },
        "lines": {
            "start": "lines:",
            "end": "links:",
        },
        "links": {
            "start": "links:",
            "end": "transformers:",
        },
        "transformers": {
            "start": "transformers:",
            "end": "# ------------------- Weather & renewables",
        },
        "atlite": {
            "start": "atlite:",
            "end": "renewable:",
        },
        "renewable_onwind": {
            "start": "onwind:",
            "end": "offwind-ac:",
        },
        "renewable_offwind-ac": {
            "start": "offwind-ac:",
            "end": "offwind-dc:",
        },
        "renewable_offwind-dc": {
            "start": "offwind-dc:",
            "end": "solar:",
        },
        "renewable_solar": {
            "start": "solar:",
            "end": "hydro:",
        },
        "renewable_hydro": {
            "start": "hydro:",
            "end": "csp:",
        },
        "renewable_csp": {
            "start": "csp:",
            "end": "# ------------------- Costs & emissions",
        },
        "costs": {
            "start": "costs:",
            "end": "co2:",
        },
        "co2": {
            "start": "co2:",
            "end": "# ------------------- Uncertainty",
        },
        "monte_carlo": {
            "start": "monte_carlo:",
            "end": "# ------------------- SECTOR OPTIONS",
        },
        "solving_options": ["solving", "options"],
        "solving_solver": ["solving", "solver"],
        "plotting": ["plotting"],
        "policy_config": ["policy_config"],
        "demand_data": ["demand_data"],
        "export": ["export"],
        "sector_toplevel": {
            "start": "sector:",
            "end": "# ------------------- HEAT SECTOR",
        },
        "sector_heat": {
            "start": "# ------------------- HEAT SECTOR",
            "end": "# ------------------- LAND TRANSPORT SECTOR",
        },
        "sector_land_transport": {
            "start": "# ------------------- LAND TRANSPORT SECTOR",
            "end": "# ------------------- BIOMASS SECTOR",
        },
        "sector_biomass": {
            "start": "# ------------------- BIOMASS SECTOR",
            "end": "# ------------------- ELECTRICITY DISTRIBUTION GRID",
        },
        "sector_electricity_distribution_grid": {
            "start": "# ------------------- ELECTRICITY DISTRIBUTION GRID",
            "end": "# ------------------- SHIPPING & AVIATION SECTOR",
        },
        "sector_shipping_aviation": {
            "start": "# ------------------- SHIPPING & AVIATION SECTOR",
            "end": "# ------------------- CCUS & CONVERSION OPTIONS",
        },
        "sector_ccus": {
            "start": "# ------------------- CCUS & CONVERSION OPTIONS",
            "end": "# ------------------- INDUSTRY OPTIONS",
        },
        "sector_industry": {
            "start": "# ------------------- INDUSTRY OPTIONS",
            "end": "# ------------------- POWERPLANTS OPTIONS",
        },
        "sector_powerplants": {
            "start": "# ------------------- POWERPLANTS OPTIONS",
            "end": "solving:",
        },
    }

    for snippet_name, keys in sections.items():
        if "start" in keys and "end" in keys:
            # Special handling - extract from "start" up to "end"
            lines = config_content.split("\n")
            start_idx = None
            end_idx = None
            for i, line in enumerate(lines):
                if line.strip().startswith(keys["start"]):
                    start_idx = i
                if start_idx is not None and line.strip().startswith(keys["end"]):
                    end_idx = i
                    break

            if start_idx is not None and end_idx is not None:
                section_text = "\n".join(lines[start_idx:end_idx])
            else:
                section_text = ""
        elif len(keys) == 1:
            # Regular top-level section
            section_text = extract_yaml_section(config_content, keys[0])
        else:
            # Nested section (e.g., renewable/onwind)
            parent_section = extract_yaml_section(config_content, keys[0])
            section_text = extract_yaml_section(parent_section, keys[1])

        if section_text.strip():
            output_file = snippets_dir / f"{snippet_name}.yaml"
            with open(output_file, "w") as f:
                # Ensure trailing newline for end-of-file-fixer pre-commit hook
                f.write(section_text.rstrip() + "\n")
            print(f"Created {output_file}")
        else:
            print(f"Warning: No content found for {snippet_name}")


if __name__ == "__main__":
    main()
