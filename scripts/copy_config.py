# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
The script copies the workflow configuration and selected core workflow files
to an output directory.

This rule is primarily used for archiving and reproducibility purposes. It
stores the configuration file and key workflow scripts alongside model outputs,
allowing the exact settings and implementation used for a simulation run to be
tracked and reproduced at a later stage.

Outputs
-------

The following files are copied:

- ``{folder}/config.yaml``: copy of the workflow configuration
- ``{folder}/Snakefile``: copy of the workflow definition
- ``{folder}/solve_network.py``: copy of the network solving script
- ``{folder}/prepare_sector_network.py``: copy of the sector network preparation script

"""
import os
from shutil import copy

from _helpers import BASE_DIR

files_to_copy = {
    os.path.join(BASE_DIR, "./config.yaml"): "config.yaml",
    os.path.join(BASE_DIR, "./Snakefile"): "Snakefile",
    os.path.join(BASE_DIR, "./scripts/solve_network.py"): "solve_network.py",
    os.path.join(
        BASE_DIR, "./scripts/prepare_sector_network.py"
    ): "prepare_sector_network.py",
}

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("copy_config")

    directory = snakemake.output["folder"]
    for f, name in files_to_copy.items():
        copy(f, directory + "/" + name)
