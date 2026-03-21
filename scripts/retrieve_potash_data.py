# SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Retrieve the USGS Potash GIS dataset.

This script downloads the Potash GIS dataset from the USGS repository
and extracts it to the local data directory. The dataset contains
global potash tract geometries which are used to estimate underground
salt cavern hydrogen storage potentials.

Source:
https://pubs.usgs.gov/sir/2010/5090/
"""

import os
import zipfile
from pathlib import Path

import requests
import shutil


def download_file(url, dest):
    response = requests.get(url)
    response.raise_for_status()
    with open(dest, "wb") as f:
        f.write(response.content)


if __name__ == "__main__":
    output_shp = Path(snakemake.output.shp)
    download_dir = output_shp.parent.parent
    zip_path = download_dir / "PotashGIS.zip"

    os.makedirs(download_dir, exist_ok=True)

    # download
    download_file(snakemake.params.url, zip_path)

    # extract
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(download_dir)

    zip_path.unlink()

    shp_file = None
    for root, _, files in os.walk(download_dir):
        for file in files:
            if file == "PotashTracts.shp":
                shp_file = Path(root) / file
                break
        if shp_file:
            break

    if shp_file is None:
        raise FileNotFoundError("PotashTracts.shp not found after extraction.")

    # ensure output dir exists
    output_shp.parent.mkdir(parents=True, exist_ok=True)

    # copy ALL shapefile components
    for ext in [".shp", ".shx", ".dbf", ".prj", ".cpg"]:
        src = shp_file.with_suffix(ext)
        if src.exists():
            shutil.copy(src, output_shp.with_suffix(ext))

    if not output_shp.exists():
        raise FileNotFoundError(f"Failed to create output: {output_shp}")


