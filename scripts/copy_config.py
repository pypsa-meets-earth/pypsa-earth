# -*- coding: utf-8 -*-
import os
from shutil import copy

files_to_copy = {
    "./config.yaml": "config.yaml",
    "./Snakefile": "Snakefile",
    "./scripts/solve_network.py": "solve_network.py",
    "./scripts/prepare_sector_network.py": "prepare_sector_network.py",
    "./config.pypsa-earth.yaml": "config.pypsa-earth.yaml",
}

if __name__ == "__main__":
    if "snakemake" not in globals():
        from helpers import mock_snakemake

        snakemake = mock_snakemake("copy_config")

    directory = snakemake.output["folder"]
    for f, name in files_to_copy.items():
        copy(f, directory + "/" + name)
