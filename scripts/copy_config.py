# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

from shutil import copy

from _helpers import BASE_DIR, get_path, mock_snakemake

files_to_copy = {
    get_path(BASE_DIR, "./config.yaml"): "config.yaml",
    get_path(BASE_DIR, "./Snakefile"): "Snakefile",
    get_path(BASE_DIR, "./scripts/solve_network.py"): "solve_network.py",
    get_path(BASE_DIR, "./scripts/prepare_sector_network.py"): "prepare_sector_network.py",
}

if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake("copy_config")

    directory = snakemake.output["folder"]
    for f, name in files_to_copy.items():
        copy(f, directory + "/" + name)
