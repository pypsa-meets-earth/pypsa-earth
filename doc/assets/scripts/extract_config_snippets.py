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
            # A non-empty line at the same or lower indentation starts the
            # next sibling section. This includes separator comments, which
            # belong to the following section rather than the current snippet.
            if line.strip():
                current_indent = len(line) - len(line.lstrip())
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
        "data_retrieval": ["enable"],
        "run": ["run"],
        "scenario": ["scenario"],
        "snapshots": ["snapshots"],
        "crs": ["crs"],
        "augmented_line_connection": ["augmented_line_connection"],
        "cluster_options": ["cluster_options"],
        "build_shape_options": ["build_shape_options"],
        "subregion": ["subregion"],
        "natura": ["natura"],
        "osm": ["osm"],
        "base_network": ["base_network"],
        "load_options": ["load_options"],
        "electricity": ["electricity"],
        "lines": ["lines"],
        "links": ["links"],
        "transformers": ["transformers"],
        "atlite": ["atlite"],
        "renewable_onwind": ["renewable", "onwind"],
        "renewable_offwind-ac": ["renewable", "offwind-ac"],
        "renewable_offwind-dc": ["renewable", "offwind-dc"],
        "renewable_solar": ["renewable", "solar"],
        "renewable_hydro": ["renewable", "hydro"],
        "renewable_csp": ["renewable", "csp"],
        "costs": ["costs"],
        "co2": ["co2"],
        "monte_carlo": ["monte_carlo"],
        "solving_solver": ["solving", "solver"],
        "solving_options": {
            "start": "# ------------------- Optimization options",
            "end": "# ------------------- Solver presets",
        },
        "plotting": ["plotting"],
        "policy_config": ["policy_config"],
        "export": ["export"],
        "demand_data": ["demand_data"],
        "custom_data": ["custom_data"],
        "existing_capacities": ["existing_capacities"],
        "sector_toplevel": {
            "start": "sector:",
            "end": "# ------------------- Heat sector",
        },
        "sector_heat": {
            "start": "# ------------------- Heat sector",
            "end": "# ------------------- Land transport sector",
        },
        "sector_land_transport": {
            "start": "# ------------------- Land transport sector",
            "end": "# ------------------- Biomass sector",
        },
        "sector_biomass": {
            "start": "# ------------------- Biomass sector",
            "end": "# ------------------- Electricity distribution grid",
        },
        "sector_electricity_distribution_grid": {
            "start": "# ------------------- Electricity distribution grid",
            "end": "# ------------------- Shipping & aviation sector",
        },
        "sector_shipping_aviation": {
            "start": "# ------------------- Shipping & aviation sector",
            "end": "# ------------------- CCUS & conversion options",
        },
        "sector_ccus": {
            "start": "# ------------------- CCUS & conversion options",
            "end": "# ------------------- Industry options",
        },
        "sector_industry": {
            "start": "# ------------------- Industry options",
            "end": "# ------------------- Powerplants options",
        },
        "sector_powerplants": {
            "start": "# ------------------- Powerplants options",
            "end": "# =================== SOLVING",
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
