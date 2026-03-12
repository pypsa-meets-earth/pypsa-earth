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


# # load tutorial hydrobasins bundle for Africa only
# bundle_tutorial_hydrobasins:
#   countries: [Africa]
#   tutorial: true
#   category: hydrobasins
#   destination: "data/hydrobasins"
#   urls:
#     hydrobasins:
#       base_url: https://data.hydrosheds.org/file/HydroBASINS/standard/
#       suffixes: ["af"]
#   unzip: true
#   output:
#   - data/hydrobasins/hybas_world.shp

# # global data for hydrobasins
# bundle_hydrobasins:
#   countries: [Earth]
#   tutorial: false
#   category: hydrobasins
#   destination: "data/hydrobasins"
#   urls:
#     hydrobasins:
#       base_url: https://data.hydrosheds.org/file/HydroBASINS/standard/
#       suffixes: ["af", "ar", "as", "au", "eu", "gr", "na", "sa", "si"]
#   unzip: true
#   output:
#   - data/hydrobasins/hybas_world.shp


if (HYDROBASINS_DATASET := dataset_version("hydrobasins"))["source"] in ["build"]:

    """
    Function to download and unzip the data for hydrobasins from HydroBASINS database
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
            unpack_archive(input["hydro_zip"], output["unzip"])

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
            import geopandas as gpd

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

            # IRENA energy statistics dataset including generation, installed capacity, heat production, and related indicators.
            # Used in build_renewable_profiles.py to derive hydropower generation potentials.
            # Original source: https://www.irena.org/-/media/Files/IRENA/Agency/Publication/2025/Jul/IRENA_Statistics_Extract_2025H2.xlsx




#   bundle_irena_statistics:
#     countries: [Earth]
#     tutorial: false
#     category: irena
#     destination: "data"
#     urls:
#       direct: https://github.com/pypsa-meets-earth/temporary_storage/raw/refs/heads/main/datasets/IRENA_Statistics_Extract_2025H2.xlsx
#     output:
#     - data/IRENA_Statistics_Extract_2025H2.xlsx

if (IRENA_DATASET := dataset_version("irena"))["source"] in ["primary"]:

    rule retrieve_irena_statistics:
        message:
            "Retrieving IRENA energy statistics dataset"
        input:
            irena_xlsx=storage(IRENA_DATASET["url"]),
        output:
            irena_xlsx_local=f"data/IRENA_Statistics_Extract_2025H2.xlsx",
        run:
            copy2(input["irena_xlsx"], output["irena_xlsx_local"])
