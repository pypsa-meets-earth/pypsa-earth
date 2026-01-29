#!/usr/bin/env python3
"""
Generate Snakemake DAGs/rulegraphs for the documentation.
This script runs snakemake --rulegraph and pipes the output to dot to generate SVGs.
Generated images are saved to doc/img/gen_*.svg and are gitignored.
"""

# SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import os
import subprocess
import sys
from pathlib import Path

# Configuration
# Map output text filenames to Snakemake arguments
# Key: Output filename suffix (e.g. 'rulegraph' -> gen_rulegraph.svg)
# Value: List of arguments to pass to snakemake
DAGS = {
    "rulegraph": ["--rulegraph", "solve_all_networks"],
    "rulegraph-sector": [
        "--rulegraph",
        "solve_sector_networks",
        "--configfile",
        "test/config.sector.yaml",
    ],
    "rulegraph-myopic": [
        "--rulegraph",
        "solve_sector_networks_myopic",
        "--configfile",
        "test/config.myopic.yaml",
    ],
    # Add more DAGs here if needed, e.g.:
    # "dag_solve": ["--dag", "solve_network"],
}

DOC_IMG_DIR = Path("doc/img")


def check_dependencies():
    """Check if snakemake and dot are available."""
    try:
        subprocess.run(["snakemake", "--version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: 'snakemake' not found. Please install it.")
        return False

    try:
        subprocess.run(["dot", "-V"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: 'dot' (graphviz) not found. Please install it.")
        return False

    return True


def generate_dags():
    """Generate the DAGs defined in DAGS."""
    if not DOC_IMG_DIR.exists():
        DOC_IMG_DIR.mkdir(parents=True, exist_ok=True)

    for name, args in DAGS.items():
        output_file = DOC_IMG_DIR / f"gen_{name}.svg"
        print(f"Generating {output_file}...")

        cmd_snakemake = ["snakemake", "-F", "-j", "1"] + args

        try:
            p_snakemake = subprocess.Popen(
                cmd_snakemake, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            # Read snakemake output
            snakemake_out, _ = p_snakemake.communicate()

            # Filter for DOT content (lines from "digraph" onwards)
            if "digraph" not in snakemake_out:
                print(f"Error generating {name}: No DOT graph found in output")
                print("Output was:", snakemake_out)
                continue

            dot_content = snakemake_out[snakemake_out.find("digraph") :]

            cmd_dot = ["dot", "-Tsvg", "-o", str(output_file)]
            p_dot = subprocess.Popen(
                cmd_dot,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            dot_out, dot_err = p_dot.communicate(input=dot_content)

            if p_dot.returncode != 0:
                print(f"Error generating {name}: dot failed")
                print(dot_err)
                continue

            # If output file exists and is non-empty, consider it a success
            if output_file.exists() and output_file.stat().st_size > 0:
                print(f"Successfully generated {output_file}")
            else:
                print(f"Error: {output_file} was not created or is empty.")
                if p_snakemake.returncode != 0:
                    print(f"Snakemake exit code: {p_snakemake.returncode}")
                    # Capture stderr from snakemake?
                    # We can't easily because we closed stdout? No, stderr is separate.
                    # But we didn't read it.
                    # We can try to read it if we didn't use communicate on p_snakemake.
                    # But p_snakemake stdout was piped. Stderr was PIPE.
                    # We should have read stderr in a thread or separate read.
                    # For simplicity, assume if dot failed or file missing, something went wrong.
                    print("Check console by running snakemake command manually.")

        except Exception as e:
            print(f"Failed to generate {name}: {e}")


if __name__ == "__main__":
    if check_dependencies():
        generate_dags()
    else:
        sys.exit(1)


def on_pre_build(config, **kwargs):
    """MkDocs hook to generate DAGs before the build."""
    if check_dependencies():
        generate_dags()
