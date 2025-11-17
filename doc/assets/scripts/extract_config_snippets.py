#!/usr/bin/env python3
"""
Extract YAML sections from config.default.yaml and save them as individual snippet files.
These snippets will automatically update when config.default.yaml changes.
"""

import re
from pathlib import Path

import yaml


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
        "toplevel": [
            "version",
            "tutorial",
            "logging",
            "results_dir",
            "summary_dir",
            "foresight",
            "countries",
            "enable",
            "custom_rules",
        ],
        "run": ["run"],
        "scenario": ["scenario"],
        "snapshots": ["snapshots"],
        "crs": ["crs"],
        "augmented_line_connection": ["augmented_line_connection"],
        "cluster_options": ["cluster_options"],
        "build_shape_options": ["build_shape_options"],
        "subregion": ["subregion"],
        "clean_osm_data_options": ["clean_osm_data_options"],
        "build_osm_network": ["build_osm_network"],
        "base_network": ["base_network"],
        "load_options": ["load_options"],
        "co2budget": ["co2_budget"],
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
        "monte_carlo": ["monte_carlo"],
        "solving_options": ["solving", "options"],
        "solving_solver": ["solving", "solver"],
        "plotting": ["plotting"],
    }

    for snippet_name, keys in sections.items():
        if snippet_name == "toplevel":
            # Special handling for toplevel - extract from 'version:' up to 'run:'
            lines = config_content.split("\n")
            start_idx = None
            end_idx = None
            for i, line in enumerate(lines):
                if line.strip().startswith("version:"):
                    start_idx = i
                if start_idx is not None and line.strip().startswith("run:"):
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
                f.write(section_text)
            print(f"Created {output_file}")
        else:
            print(f"Warning: No content found for {snippet_name}")


if __name__ == "__main__":
    main()
