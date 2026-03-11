# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import os
import requests
from datetime import datetime
from dateutil.relativedelta import relativedelta
from shutil import move, unpack_archive, rmtree, copy2
from zipfile import ZipFile


# Configure the default storage provider for accessing remote files using http
# and the special storage plugin for accessing Zenodo files
storage:
    provider="http",
    keep_local=True,
    retries=3,


storage cached_http:
    provider="cached-http",


if (HYDROBASINS_DATASET := dataset_version("hydrobasins"))["source"] in ["build"]:

    suffixes = ["af", "ar", "as", "au", "eu", "gr", "na", "sa", "si"]
    level = config["renewable"]["hydro"]["hydrobasins_level"]

    rule retrieve_hydrobasins:
        message:
            "Retrieving hydrobasins dataset for {wildcards.suffix}"
        input:
            hydro_zip=storage(
                f"{HYDROBASINS_DATASET['url']}" + "/hybas_{suffix}_lev01-12_v1c.zip"
            ),
        output:
            unzip=directory(
                f"{HYDROBASINS_DATASET['folder']}" + "/hybas_{suffix}_lev01-12_v1c"
            ),
            shp=f"{HYDROBASINS_DATASET['folder']}"
            + "/hybas_{suffix}_lev01-12_v1c"
            + "/hybas_{suffix}_"
            + f"lev{level:02d}_v1c.shp",
        run:
            unpack_archive(input["hydro_zip"], output["unzip"])

    rule create_hydrobasins_world:
        message:
            "Aggregate hydrobasins into single dataset"
        input:
            expand(
                f"{HYDROBASINS_DATASET['folder']}"
                + "/hybas_{suffix}_lev01-12_v1c"
                + "/hybas_{suffix}_"
                + f"lev{level:02d}_v1c.shp",
                suffix=suffixes,
            ),
        output:
            "data/hydrobasins/hybas_world.shp",
        run:
            import geopandas as gpd

            gpdf_list = []
            logger.info(f"Merging hydrobasins files into: {output}")
            for f_name in input:
                gpdf_list.append(gpd.read_file(f_name))

            merged = gpd.GeoDataFrame(pd.concat(gpdf_list)).drop_duplicates(
                subset="HYBAS_ID", ignore_index=True
            )
            merged.to_file(str(output), driver="ESRI Shapefile")
