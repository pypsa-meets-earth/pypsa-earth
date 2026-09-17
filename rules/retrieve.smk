# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

from retrieve_databundle_light import (
    datafiles_retrivedatabundle,
    get_best_bundles_in_snakemake,
)
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

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


if config["enable"].get("retrieve_databundle", True):

    bundles_to_download = get_best_bundles_in_snakemake(
        config, exclude_categories=["cutouts"]
    )

    rule retrieve_databundle_light:
        params:
            bundles_to_download=bundles_to_download,
            hydrobasins_level=config["renewable"]["hydro"]["hydrobasins_level"],
        output:  #expand(directory('{file}') if isdir('{file}') else '{file}', file=datafiles)
            expand(
                "{file}", file=datafiles_retrivedatabundle(config, bundles_to_download)
            ),
            directory("data/landcover"),
        log:
            "logs/" + RDIR + "retrieve_databundle.log",
        benchmark:
            "benchmarks/" + RDIR + "retrieve_databundle_light"
        script:
            scripts("retrieve_databundle_light.py")


if config["enable"].get("download_global_buildings", True):

    rule download_global_buildings:
        params:
            crs=config["crs"],
        output:
            "data/global_buildings/{country}_global_buildings_raw.parquet",
        script:
            scripts("download_global_buildings.py")


if config["enable"].get("download_osm_data", True):

    rule download_osm_data:
        params:
            countries=config["countries"],
        output:
            cables="resources/" + RDIR + "osm/raw/all_raw_cables.geojson",
            generators="resources/" + RDIR + "osm/raw/all_raw_generators.geojson",
            generators_csv="resources/" + RDIR + "osm/raw/all_raw_generators.csv",
            lines="resources/" + RDIR + "osm/raw/all_raw_lines.geojson",
            substations="resources/" + RDIR + "osm/raw/all_raw_substations.geojson",
        log:
            "logs/" + RDIR + "download_osm_data.log",
        benchmark:
            "benchmarks/" + RDIR + "download_osm_data"
        script:
            scripts("download_osm_data.py")


if config["enable"].get("retrieve_cutout", False):

    cutout_to_download = get_best_bundles_in_snakemake(
        config, include_categories=["cutouts"]
    )

    rule retrieve_cutout:
        params:
            bundles_to_download=cutout_to_download,
            hydrobasins_level=[],
        input:
            check=terminate_if_cutout_exists,
            onshore_shapes="resources/" + RDIR + "shapes/country_shapes.geojson",
            offshore_shapes="resources/" + RDIR + "shapes/offshore_shapes.geojson",
        output:
            "cutouts/" + CDIR + "{cutout}.nc",
        log:
            "logs/" + RDIR + "retrieve_cutout/{cutout}.log",
        benchmark:
            "benchmarks/" + RDIR + "retrieve_cutout_{cutout}"
        script:
            scripts("retrieve_databundle_light.py")


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
            costs="resources/" + RDIR + "costs_{year}.csv",
        log:
            "logs/" + RDIR + "retrieve_cost_data_{year}.log",
        resources:
            mem_mb=5000,
        run:
            move(input[0], output[0])


rule retrieve_potash_data:
    input:
        potash_zip=HTTP.remote(
            "https://pubs.usgs.gov/sir/2010/5090/s/PotashGIS.zip",
            keep_local=True,
        ),
    output:
        potash_dir=directory("data/potash_gis"),
        potash_files="data/potash_gis/PotashGIS/global_potash/Shapefiles/PotashTracts.shp",
    run:
        unpack_archive(str(input.potash_zip), output["potash_dir"])


if config["co2"]["automatic_emission"]["enable"]:

    rule retrieve_emissions:
        input:
            HTTP.remote(
                "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/v60_GHG/CO2_excl_short-cycle_org_C/v60_GHG_CO2_excl_short-cycle_org_C_1970_2018.zip",
                keep_local=True,
            ),
        output:
            edgar_folder=directory("data/co2_emissions/"),
            edgar_zip="data/co2_emissions/v60_GHG_CO2_excl_short-cycle_org_C_1970_2018.zip",
            edgar_xlsx="data/co2_emissions/v60_CO2_excl_short-cycle_org_C_1970_2018.xls",
        log:
            "logs/" + RDIR + "retrieve_emissions.log",
        run:
            move(input[0], output.edgar_zip)
            unpack_archive(output.edgar_zip, extract_dir=output.edgar_folder)


rule retrieve_us_cities_dataset:
    output:
        us_cities="data/industry/us_cities.csv",
    script:
        scripts("retrieve_us_cities_dataset.py")


rule retrieve_ammonia_dataset:
    output:
        usgs_ammonia_dataset="data/industry/USGS_ammonia_dataset.xlsx",
    script:
        scripts("retrieve_ammonia_dataset.py")
