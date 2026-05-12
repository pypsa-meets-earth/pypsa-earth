# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import os
import requests
from datetime import datetime
from dateutil.relativedelta import relativedelta
from shutil import move, unpack_archive, rmtree, copy2
from zipfile import ZipFile
from scripts._common import dataset_version


if (HYDROBASINS_DATASET := dataset_version("hydrobasins", config))["source"] in [
    "build",
    "tutorial",
]:

    """
    Rules to download and unzip the data for hydrobasins from HydroBASINS database
    available via https://www.hydrosheds.org/products/hydrobasins

    We are using data from the HydroSHEDS version 1 database
    which is © World Wildlife Fund, Inc. (2006-2022) and has been used herein under license.
    WWF has not evaluated our data pipeline and therefore gives no warranty regarding its
    accuracy, completeness, currency or suitability for any particular purpose.
    Portions of the HydroSHEDS v1 database incorporate data which are the intellectual property
    rights of © USGS (2006-2008), NASA (2000-2005), ESRI (1992-1998), CIAT (2004-2006),
    UNEP-WCMC (1993), WWF (2004), Commonwealth of Australia (2007), and Her Royal Majesty
    and the British Crown and are used under license. The HydroSHEDS v1 database and
    more information are available at https://www.hydrosheds.org.
    """
    suffixes: list = HYDROBASINS_DATASET["region"].split(" ")
    level = config["renewable"]["hydro"]["hydrobasins_level"]

    def get_hydrobasins(input, output):
        """Merge the hydrobasins files into a single shapefile

        Arguments
        ---------
        input : list
            List of input files to merge
        output : dict
            Dictionary of output files to write the merged shapefile
        """

        import geopandas as gpd
        import pandas as pd

        gpdf_list = []
        logger.info(f"Merging hydrobasins files into: {output}")
        for f_name in input:
            if f_name.endswith(".shp"):
                logger.info(f"Reading hydrobasins file: {f_name}")
                gpdf_list.append(gpd.read_file(f_name))

        merged = gpd.GeoDataFrame(pd.concat(gpdf_list)).drop_duplicates(
            subset="HYBAS_ID", ignore_index=True
        )
        merged.to_file(str(output["shp"]), driver="ESRI Shapefile")

    rule retrieve_hydrobasins:
        message:
            "Retrieving hydrobasins dataset for {wildcards.suffix}"
        input:
            hydro_zip=HTTP.remote(
                f"{HYDROBASINS_DATASET['url']}" + "/hybas_{suffix}_lev01-12_v1c.zip",
                keep_local=True,
            ),
        output:
            unzip=directory(
                f"{HYDROBASINS_DATASET['folder']}" + "/hybas_{suffix}_lev01-12_v1c"
            ),
            shp=multiext(
                f"{HYDROBASINS_DATASET['folder']}"
                + "/hybas_{suffix}_lev01-12_v1c"
                + "/hybas_{suffix}_"
                + f"lev{level:02d}_v1c",
                ".dbf",
                ".prj",
                ".sbn",
                ".sbx",
                ".shp",
                ".shp.xml",
                ".shx",
            ),
        run:
            unpack_archive(str(input["hydro_zip"]), output["unzip"])

    rule create_hydrobasins_world:
        message:
            "Aggregate hydrobasins into single dataset"
        input:
            expand(
                f"{HYDROBASINS_DATASET['folder']}"
                + "/hybas_{suffix}_lev01-12_v1c"
                + "/hybas_{suffix}_"
                + f"lev{level:02d}_v1c"
                + "{ext}",
                suffix=suffixes,
                ext=[".dbf", ".prj", ".shp", ".shx"],
            ),
        output:
            shp="data/hydrobasins/hybas_world.shp",
            other=multiext(
                "data/hydrobasins/hybas_world", ".cpg", ".dbf", ".prj", ".shx"
            ),
        run:
            get_hydrobasins(input, output)


if (IRENA_DATASET := dataset_version("irena", config))["source"] in ["primary"]:

    rule retrieve_irena_statistics:
        message:
            "Retrieving IRENA energy statistics dataset"
        input:
            irena_xlsx=HTTP.remote(IRENA_DATASET["url"], keep_local=True),
        output:
            irena_xlsx_local=f"data/IRENA_Statistics_Extract_2025H2.xlsx",
        run:
            copy2(str(input["irena_xlsx"]), output["irena_xlsx_local"])


if (LANDCOVER_DATASET := dataset_version("landcover", config))["source"] in ["primary"]:

    folder = LANDCOVER_DATASET["folder"]
    version = LANDCOVER_DATASET["version"]

    rule retrieve_landcover:
        message:
            "Retrieving landcover dataset"
        input:
            landcover_zip=HTTP.remote(LANDCOVER_DATASET["url"], keep_local=True),
        output:
            unzip=directory(f"{folder}"),
            zips=expand(
                f"{folder}/WDPA_{version}" + "_Public_shp_{index}.zip", index=[0, 1, 2]
            ),
        run:
            unpack_archive(str(input["landcover_zip"]), output["unzip"])

    rule unpack_landcover_zips:
        input:
            zip=f"{folder}/WDPA_{version}" + "_Public_shp_{index}.zip",
        output:
            dir=directory(
                f"data/landcover/world_protected_areas/WDPA_{version}"
                + "_Public_shp_{index}"
            ),
            shp=f"data/landcover/world_protected_areas/WDPA_{version}"
            + "_Public_shp_{index}/WDPA_"
            + f"{version}_Public_shp-points.shp",
        run:
            unpack_archive(str(input["zip"]), output["dir"])

    rule target_landcover:
        input:
            expand(
                f"data/landcover/world_protected_areas/WDPA_{version}_"
                + "Public_shp_{index}"
                + f"/WDPA_{version}_Public_shp-points.shp",
                index=[0, 1, 2],
            ),
