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


country_data = config["costs"].get("country_specific_data", "")
countries = config.get("countries", [])

if country_data and countries == [country_data]:
    cost_directory = f"{country_data}/"
elif country_data:
    cost_directory = f"{country_data}/"
    warnings.warn(
        f"'country_specific_data' is set to '{country_data}', but 'countries' is {countries}. Make sure the '{country_data}' directory exists and that this is intentional."
    )
else:
    cost_directory = ""


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


rule download_custom_powerplants:
    input:
        url=HTTP.remote(
            "https://sandbox.zenodo.org/records/499641/files/custom_powerplants.csv",
            keep_local=True,
            additional_request_string="?download=1",
        ),
    output:
        "data/custom_powerplants.csv",
    log:
        "logs/download_custom_powerplants.log",
    run:
        copyfile(str(input["url"]), output[0])


rule download_interconnection_data:
    input:
        substations=HTTP.remote(
            "https://sandbox.zenodo.org/records/471583/files/zm_substations.csv",
            keep_local=True,
            additional_request_string="?download=1",
        ),
        links=HTTP.remote(
            "https://sandbox.zenodo.org/records/471583/files/sapp_links.csv",
            keep_local=True,
            additional_request_string="?download=1",
        ),
        countries=HTTP.remote(
            "https://sandbox.zenodo.org/records/471583/files/sapp_countries.csv",
            keep_local=True,
            additional_request_string="?download=1",
        ),
    output:
        substations="data/zm_substations.csv",
        links="data/sapp_links.csv",
        countries="data/sapp_countries.csv",
    log:
        "logs/download_interconnection_data.log",
    run:
        copyfile(str(input["substations"]), output["substations"])
        copyfile(str(input["links"]), output["links"])
        copyfile(str(input["countries"]), output["countries"])


rule download_line_types:
    input:
        url=HTTP.remote(
            "https://sandbox.zenodo.org/records/473405/files/pypsa_line_types%20%281%29.csv",
            keep_local=True,
        ),
    output:
        "data/line_types.csv",
    log:
        "logs/download_line_types.log",
    run:
        copyfile(str(input["url"]), output[0])


rule retrieve_mining_data:
    input:
        provincial_demand=HTTP.remote(
            "https://sandbox.zenodo.org/records/495635/files/zambia_provincial_mining_demand.csv",
            keep_local=True,
            additional_request_string="?download=1",
        ),
        mining_polygons=HTTP.remote(
            "https://sandbox.zenodo.org/records/495635/files/zambia_pangaea_mining_polygons.csv",
            keep_local=True,
            additional_request_string="?download=1",
        ),
    output:
        provincial_demand="data/mining/zambia_provincial_mining_demand.csv",
        mining_polygons="data/mining/zambia_pangaea_mining_polygons.csv",
    log:
        "logs/retrieve_mining_data.log",
    run:
        import os

        os.makedirs("data/mining", exist_ok=True)
        copyfile(str(input["provincial_demand"]), output["provincial_demand"])
        copyfile(str(input["mining_polygons"]), output["mining_polygons"])


if config["enable"].get("retrieve_cost_data", True):

    rule retrieve_cost_data:
        params:
            version=config["costs"]["technology_data_version"],
        input:
            HTTP.remote(
                f"raw.githubusercontent.com/PyPSA/technology-data/{config['costs']['technology_data_version']}/outputs/{cost_directory}"
                + "costs_{year}.csv",
                keep_local=True,
            ),
        output:
            "resources/" + RDIR + "costs_{year}.csv",
        log:
            "logs/" + RDIR + "retrieve_cost_data_{year}.log",
        resources:
            mem_mb=5000,
        run:
            move(input[0], output[0])


if (HYDRO_PROFILE_DATASET := dataset_version("hydro_profile", config))["source"] in [
    "primary",
    "tutorial",
]:

    region = HYDRO_PROFILE_DATASET["region"]
    source = HYDRO_PROFILE_DATASET["source"]

    rule retrieve_hydro_profile:
        message:
            "Retrieving hydro profile dataset for {region} and {source}"
        input:
            hydro_profile_nc=HTTP.remote(
                HYDRO_PROFILE_DATASET["url"],
                keep_local=True,
                additional_request_string="?download=1",
            ),
        output:
            f"data/hydro_profiles/glofas_profile.nc",
        run:
            copy2(str(input[0]), output[0])


if (NATURA_EARTH_DATASET := dataset_version("natura_earth", config))["source"] in [
    "primary",
    "tutorial",
    "archive",
]:

    source = NATURA_EARTH_DATASET["source"]

    rule retrieve_natura_earth:
        message:
            "Retrieving Natura Earth dataset for {source}"
        input:
            natura_zip=HTTP.remote(
                NATURA_EARTH_DATASET["url"],
                keep_local=True,
                additional_request_string="?download=1",
            ),
        output:
            unzip=directory(f"data/natura_earth/{source}"),
            tiff=f"data/natura_earth/{source}" + "/natura.tiff",
            shp=f"data/natura/natura.tiff",
        run:
            unpack_archive(str(input["natura_zip"]), output["unzip"])
            copy2(os.path.join(output["tiff"]), output["shp"])


if (ERA5_CUTOUT := dataset_version("cutout-era5", config))["source"] in [
    "primary",
    "tutorial",
]:
    year = int(float(ERA5_CUTOUT["year"]))
    region = ERA5_CUTOUT["region"]

    rule retrieve_era5_cutout:
        message:
            f"Retrieving ERA5 cutout for region {region} ({year})"
        input:
            cutout=HTTP.remote(
                ERA5_CUTOUT["url"],
                keep_local=True,
                additional_request_string="?download=1",
            ),
        output:
            f"cutouts/{CDIR}cutout-{year}-era5.nc",
        log:
            f"logs/{RDIR}retrieve_era5_cutout.log",
        benchmark:
            f"benchmarks/{RDIR}retrieve_era5_cutout"
        run:
            copy2(str(input[0]), output[0])
