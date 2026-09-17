# SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Generate cutout coverage plots before MkDocs collects documentation files."""

import subprocess
import sys
from pathlib import Path


def on_pre_build(config):
    """Create documentation plots from the checked-in cutout metadata CSV."""

    repository = Path(__file__).resolve().parents[3]
    subprocess.run(
        [
            sys.executable,
            repository / "scripts/non_workflow/extract_cutout_bundle_info.py",
            "--skip-csv",
            "--csv-path",
            repository / "doc/user-guide/data/cutout_bundle_info.csv",
            "--output-dir",
            repository / "doc/user-guide/images",
        ],
        cwd=repository,
        check=True,
    )
