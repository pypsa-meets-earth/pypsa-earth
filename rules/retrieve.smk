# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

from retrieve_databundle_light import (
    datafiles_retrivedatabundle,
    get_best_bundles_in_snakemake,
    get_databundle_categories,
)

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

    # Collect all databundle categories, except cutouts that are treated separately
    databundle_categories = get_databundle_categories(
        config["databundles"], exclude_categories=["cutouts"]
    )

    # get the best bundles for each category
    bundle_dict = {
        category: get_best_bundles_in_snakemake(config, include_categories=[category])
        for category in databundle_categories
    }

    # List selected bundles for each category
    print("\n====================================================")
    print("Selected bundles for each category:")
    for category, bundles in bundle_dict.items():
        if len(bundles) >= 1:
            print(f"\t{category}: {bundles}")
    print("====================================================\n")

    for category, bundles in bundle_dict.items():
        output_files = datafiles_retrivedatabundle(config, bundles)
        contains_landcover = any([f for f in output_files if "landcover" in f])

        if not bundles:
            continue

        rule:
            name:
                f"retrieve_databundle_{category}"
            params:
                tutorial=config["tutorial"],
                bundles_to_download=bundles,
                hydrobasins_level=config["renewable"]["hydro"]["hydrobasins_level"],
            output:
                expand("{file}", file=output_files),
                branch(contains_landcover, directory("data/landcover")),
            log:
                "logs/" + RDIR + f"retrieve_databundle_{category}.log",
            benchmark:
                "benchmarks/" + RDIR + f"retrieve_databundle_{category}"
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

    print("Cutout to download: ", cutout_to_download)

    rule retrieve_cutout:
        params:
            tutorial=config["tutorial"],
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
            storage.http(
                f"https://raw.githubusercontent.com/PyPSA/technology-data/{config['costs']['technology_data_version']}/outputs/{cost_directory}"
                + "costs_{planning_horizons}.csv",
                keep_local=True,
            ),
        output:
            costs="resources/" + RDIR + "costs_{planning_horizons}.csv",
        log:
            "logs/" + RDIR + "retrieve_cost_data_{planning_horizons}.log",
        resources:
            mem_mb=5000,
        run:
            move(input[0], output[0])


rule retrieve_potash_data:
    input:
        potash_zip=storage.http(
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
            storage.http(
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


rule retrieve_urban_percent:
    params:
        url_urban_percent="https://zenodo.org/records/23242372/files/US_PopTotal.csv.7z?download=1",
    output:
        urban_percent_raw="data/urban_percent_raw.csv",
    script:
        scripts("retrieve_urban_percent.py")
