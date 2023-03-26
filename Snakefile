# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later

import sys

sys.path.append("./scripts")

from os.path import normpath, exists, isdir
from shutil import copyfile, move

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

from scripts._helpers import create_country_list
from scripts.build_demand_profiles import get_load_paths_gegis
from scripts.retrieve_databundle_light import datafiles_retrivedatabundle
from pathlib import Path

HTTP = HTTPRemoteProvider()

if "config" not in globals() or not config:  # skip when used as sub-workflow
    if not exists("config.yaml"):
        copyfile("config.tutorial.yaml", "config.yaml")

    configfile: "config.yaml"


configfile: "configs/bundle_config.yaml"


DEFAULT_CONFIG = "config.tutorial.yaml" if config["tutorial"] else "config.default.yaml"


# convert country list according to the desired region
config["countries"] = create_country_list(config["countries"])


# create a list of iteration steps, required to solve the experimental design
# each value is used as wildcard input e.g. solution_{unc}
config["scenario"]["unc"] = [
    f"m{i}" for i in range(config["monte_carlo"]["options"]["samples"])
]

run = config.get("run", {})
RDIR = run["name"] + "/" if run.get("name") else ""
CDIR = RDIR if not run.get("shared_cutouts") else ""

load_data_paths = get_load_paths_gegis("data", config)
if config["enable"].get("retrieve_cost_data", True):
    COSTS = "resources/" + RDIR + "costs.csv"
else:
    COSTS = "data/costs.csv"
ATLITE_NPROCESSES = config["atlite"].get("nprocesses", 10)


wildcard_constraints:
    simpl="[a-zA-Z0-9]*|all",
    clusters="[0-9]+m?|all",
    ll="(v|c)([0-9\.]+|opt|all)|all",
    opts="[-+a-zA-Z0-9\.]*",
    unc="[-+a-zA-Z0-9\.]*",


rule clean:
    shell:
        "snakemake -j 1 solve_all_networks --delete-all-output"


rule run_tests:
    output:
        touch("tests.done"),
    run:
        import os

        shell("snakemake --cores all build_test_configs")
        directory = "test/tmp"  # assign directory
        for filename in os.scandir(directory):  # iterate over files in that directory
            if filename.is_file():
                print(filename.path)
                shell("cp {filename.path} config.yaml")
                if "monte" in filename.name:
                    shell("snakemake --cores all solve_all_networks_monte --forceall")
                else:
                    shell("snakemake --cores all solve_all_networks --forceall")
        print("Tests are successful.")


rule solve_all_networks:
    input:
        expand(
            "results/" + RDIR + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            **config["scenario"]
        ),


rule plot_all_p_nom:
    input:
        expand(
            "results/"
            + RDIR
            + "plots/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_p_nom.{ext}",
            **config["scenario"],
            ext=["png", "pdf"]
        ),


rule make_all_summaries:
    input:
        expand(
            "results/"
            + RDIR
            + "summaries/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}",
            **config["scenario"],
            country=["all"] + config["countries"]
        ),


rule plot_all_summaries:
    input:
        expand(
            "results/"
            + RDIR
            + "plots/summary_{summary}_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}.{ext}",
            summary=["energy", "costs"],
            **config["scenario"],
            country=["all"] + config["countries"],
            ext=["png", "pdf"]
        ),


if config["enable"].get("retrieve_databundle", True):

    rule retrieve_databundle_light:
        output:  #expand(directory('{file}') if isdir('{file}') else '{file}', file=datafiles)
            expand("{file}", file=datafiles_retrivedatabundle(config)),
            directory("data/landcover"),
        log:
            "logs/" + RDIR + "retrieve_databundle.log",
        benchmark:
            "benchmarks/" + RDIR + "retrieve_databundle_light"
        script:
            "scripts/retrieve_databundle_light.py"


if config["enable"].get("download_osm_data", True):

    rule download_osm_data:
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
            "scripts/download_osm_data.py"


rule clean_osm_data:
    input:
        cables="resources/" + RDIR + "osm/raw/all_raw_cables.geojson",
        generators="resources/" + RDIR + "osm/raw/all_raw_generators.geojson",
        lines="resources/" + RDIR + "osm/raw/all_raw_lines.geojson",
        substations="resources/" + RDIR + "osm/raw/all_raw_substations.geojson",
        country_shapes="resources/" + RDIR + "shapes/country_shapes.geojson",
        offshore_shapes="resources/" + RDIR + "shapes/offshore_shapes.geojson",
        africa_shape="resources/" + RDIR + "shapes/africa_shape.geojson",
    output:
        generators="resources/" + RDIR + "osm/clean/all_clean_generators.geojson",
        generators_csv="resources/" + RDIR + "osm/clean/all_clean_generators.csv",
        lines="resources/" + RDIR + "osm/clean/all_clean_lines.geojson",
        substations="resources/" + RDIR + "osm/clean/all_clean_substations.geojson",
    log:
        "logs/" + RDIR + "clean_osm_data.log",
    benchmark:
        "benchmarks/" + RDIR + "clean_osm_data"
    script:
        "scripts/clean_osm_data.py"


rule build_osm_network:
    input:
        generators="resources/" + RDIR + "osm/clean/all_clean_generators.geojson",
        lines="resources/" + RDIR + "osm/clean/all_clean_lines.geojson",
        substations="resources/" + RDIR + "osm/clean/all_clean_substations.geojson",
        country_shapes="resources/" + RDIR + "shapes/country_shapes.geojson",
    output:
        lines="resources/" + RDIR + "base_network/all_lines_build_network.csv",
        converters="resources/" + RDIR + "base_network/all_converters_build_network.csv",
        transformers="resources/"
        + RDIR
        + "base_network/all_transformers_build_network.csv",
        substations="resources/" + RDIR + "base_network/all_buses_build_network.csv",
    log:
        "logs/" + RDIR + "build_osm_network.log",
    benchmark:
        "benchmarks/" + RDIR + "build_osm_network"
    script:
        "scripts/build_osm_network.py"


rule build_shapes:
    input:
        # naturalearth='data/bundle/naturalearth/ne_10m_admin_0_countries.shp',
        # eez='data/bundle/eez/World_EEZ_v8_2014.shp',
        # nuts3='data/bundle/NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp',
        # nuts3pop='data/bundle/nama_10r_3popgdp.tsv.gz',
        # nuts3gdp='data/bundle/nama_10r_3gdp.tsv.gz',
        eez="data/eez/eez_v11.gpkg",
    output:
        country_shapes="resources/" + RDIR + "shapes/country_shapes.geojson",
        offshore_shapes="resources/" + RDIR + "shapes/offshore_shapes.geojson",
        africa_shape="resources/" + RDIR + "shapes/africa_shape.geojson",
        gadm_shapes="resources/" + RDIR + "shapes/gadm_shapes.geojson",
    log:
        "logs/" + RDIR + "build_shapes.log",
    benchmark:
        "benchmarks/" + RDIR + "build_shapes"
    threads: 1
    resources:
        mem_mb=500,
    script:
        "scripts/build_shapes.py"


rule base_network:
    input:
        osm_buses="resources/" + RDIR + "base_network/all_buses_build_network.csv",
        osm_lines="resources/" + RDIR + "base_network/all_lines_build_network.csv",
        osm_converters="resources/"
        + RDIR
        + "base_network/all_converters_build_network.csv",
        osm_transformers="resources/"
        + RDIR
        + "base_network/all_transformers_build_network.csv",
        country_shapes="resources/" + RDIR + "shapes/country_shapes.geojson",
        offshore_shapes="resources/" + RDIR + "shapes/offshore_shapes.geojson",
    output:
        "networks/" + RDIR + "base.nc",
    log:
        "logs/" + RDIR + "base_network.log",
    benchmark:
        "benchmarks/" + RDIR + "base_network"
    threads: 1
    resources:
        mem_mb=500,
    script:
        "scripts/base_network.py"


rule build_bus_regions:
    input:
        country_shapes="resources/" + RDIR + "shapes/country_shapes.geojson",
        offshore_shapes="resources/" + RDIR + "shapes/offshore_shapes.geojson",
        base_network="networks/" + RDIR + "base.nc",
        #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson",
        #using this line instead of the following will test updated gadm shapes for MA.
        #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
        #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
        gadm_shapes="resources/" + RDIR + "shapes/gadm_shapes.geojson",
    output:
        regions_onshore="resources/" + RDIR + "bus_regions/regions_onshore.geojson",
        regions_offshore="resources/" + RDIR + "bus_regions/regions_offshore.geojson",
    log:
        "logs/" + RDIR + "build_bus_regions.log",
    benchmark:
        "benchmarks/" + RDIR + "build_bus_regions"
    threads: 1
    resources:
        mem_mb=1000,
    script:
        "scripts/build_bus_regions.py"


if config["enable"].get("build_cutout", False):

    rule build_cutout:
        input:
            onshore_shapes="resources/" + RDIR + "shapes/country_shapes.geojson",
            offshore_shapes="resources/" + RDIR + "shapes/offshore_shapes.geojson",
        output:
            "cutouts/" + CDIR + "{cutout}.nc",
        log:
            "logs/" + RDIR + "build_cutout/{cutout}.log",
        benchmark:
            "benchmarks/" + RDIR + "build_cutout_{cutout}"
        threads: ATLITE_NPROCESSES
        resources:
            mem_mb=ATLITE_NPROCESSES * 1000,
        script:
            "scripts/build_cutout.py"


if config["enable"].get("build_natura_raster", False):

    rule build_natura_raster:
        input:
            shapefiles_land="data/landcover",
            cutouts=expand("cutouts/" + CDIR + "{cutouts}.nc", **config["atlite"]),
        output:
            "resources/" + RDIR + "natura.tiff",
        log:
            "logs/" + RDIR + "build_natura_raster.log",
        benchmark:
            "benchmarks/" + RDIR + "build_natura_raster"
        script:
            "scripts/build_natura_raster.py"


if not config["enable"].get("build_natura_raster", False):

    rule copy_defaultnatura_tiff:
        input:
            "data/natura.tiff",
        output:
            "resources/" + RDIR + "natura.tiff",
        run:
            import shutil

            shutil.copyfile(input[0], output[0])


if config["enable"].get("retrieve_cost_data", True):

    rule retrieve_cost_data:
        input:
            HTTP.remote(
                f"raw.githubusercontent.com/PyPSA/technology-data/{config['costs']['version']}/outputs/costs_{config['costs']['year']}.csv",
                keep_local=True,
            ),
        output:
            COSTS,
        log:
            "logs/" + RDIR + "retrieve_cost_data.log",
        resources:
            mem_mb=5000,
        run:
            move(input[0], output[0])


rule build_demand_profiles:
    input:
        base_network="networks/" + RDIR + "base.nc",
        regions="resources/" + RDIR + "bus_regions/regions_onshore.geojson",
        load=load_data_paths,
        #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson", 
        #using this line instead of the following will test updated gadm shapes for MA.
        #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
        #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
        gadm_shapes="resources/" + RDIR + "shapes/gadm_shapes.geojson",
    output:
        "resources/" + RDIR + "demand_profiles.csv",
    log:
        "logs/" + RDIR + "build_demand_profiles.log",
    benchmark:
        "benchmarks/" + RDIR + "build_demand_profiles"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/build_demand_profiles.py"


rule build_renewable_profiles:
    input:
        natura="resources/" + RDIR + "natura.tiff",
        copernicus="data/copernicus/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        gebco="data/gebco/GEBCO_2021_TID.nc",
        country_shapes="resources/" + RDIR + "shapes/country_shapes.geojson",
        offshore_shapes="resources/" + RDIR + "shapes/offshore_shapes.geojson",
        hydro_capacities="data/hydro_capacities.csv",
        eia_hydro_generation="data/eia_hydro_annual_generation.csv",
        powerplants="resources/" + RDIR + "powerplants.csv",
        regions=lambda w: (
            "resources/" + RDIR + "bus_regions/regions_onshore.geojson"
            if w.technology in ("onwind", "solar", "hydro")
            else "resources/" + RDIR + "bus_regions/regions_offshore.geojson"
        ),
        cutout=lambda w: "cutouts/"
        + CDIR
        + config["renewable"][w.technology]["cutout"]
        + ".nc",
    output:
        profile="resources/" + RDIR + "renewable_profiles/profile_{technology}.nc",
    log:
        "logs/" + RDIR + "build_renewable_profile_{technology}.log",
    benchmark:
        "benchmarks/" + RDIR + "build_renewable_profiles_{technology}"
    threads: ATLITE_NPROCESSES
    resources:
        mem_mb=ATLITE_NPROCESSES * 5000,
    script:
        "scripts/build_renewable_profiles.py"


rule build_powerplants:
    input:
        base_network="networks/" + RDIR + "base.nc",
        pm_config="configs/powerplantmatching_config.yaml",
        custom_powerplants="data/custom_powerplants.csv",
        osm_powerplants="resources/" + RDIR + "osm/clean/all_clean_generators.csv",
        #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson",
        #using this line instead of the following will test updated gadm shapes for MA.
        #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
        #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
        gadm_shapes="resources/" + RDIR + "shapes/gadm_shapes.geojson",
    output:
        powerplants="resources/" + RDIR + "powerplants.csv",
        powerplants_osm2pm="resources/" + RDIR + "powerplants_osm2pm.csv",
    log:
        "logs/" + RDIR + "build_powerplants.log",
    benchmark:
        "benchmarks/" + RDIR + "build_powerplants"
    threads: 1
    resources:
        mem_mb=500,
    script:
        "scripts/build_powerplants.py"


rule add_electricity:
    input:
        **{
            f"profile_{tech}": "resources/"
            + RDIR
            + f"renewable_profiles/profile_{tech}.nc"
            for tech in config["renewable"]
            if tech in config["electricity"]["renewable_carriers"]
        },
        **{
            f"conventional_{carrier}_{attr}": fn
            for carrier, d in config.get("conventional", {None: {}}).items()
            for attr, fn in d.items()
            if str(fn).startswith("data/")
        },
        base_network="networks/" + RDIR + "base.nc",
        tech_costs=COSTS,
        powerplants="resources/" + RDIR + "powerplants.csv",
        #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson", 
        #using this line instead of the following will test updated gadm shapes for MA.
        #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
        #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
        gadm_shapes="resources/" + RDIR + "shapes/gadm_shapes.geojson",
        hydro_capacities="data/hydro_capacities.csv",
        demand_profiles="resources/" + RDIR + "demand_profiles.csv",
    output:
        "networks/" + RDIR + "elec.nc",
    log:
        "logs/" + RDIR + "add_electricity.log",
    benchmark:
        "benchmarks/" + RDIR + "add_electricity"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/add_electricity.py"


rule simplify_network:
    input:
        network="networks/" + RDIR + "elec.nc",
        tech_costs=COSTS,
        regions_onshore="resources/" + RDIR + "bus_regions/regions_onshore.geojson",
        regions_offshore="resources/" + RDIR + "bus_regions/regions_offshore.geojson",
    output:
        network="networks/" + RDIR + "elec_s{simpl}.nc",
        regions_onshore="resources/"
        + RDIR
        + "bus_regions/regions_onshore_elec_s{simpl}.geojson",
        regions_offshore="resources/"
        + RDIR
        + "bus_regions/regions_offshore_elec_s{simpl}.geojson",
        busmap="resources/" + RDIR + "bus_regions/busmap_elec_s{simpl}.csv",
        connection_costs="resources/"
        + RDIR
        + "bus_regions/connection_costs_s{simpl}.csv",
    log:
        "logs/" + RDIR + "simplify_network/elec_s{simpl}.log",
    benchmark:
        "benchmarks/" + RDIR + "simplify_network/elec_s{simpl}"
    threads: 1
    resources:
        mem_mb=4000,
    script:
        "scripts/simplify_network.py"


if config["augmented_line_connection"].get("add_to_snakefile", False) == True:

    rule cluster_network:
        input:
            network="networks/" + RDIR + "elec_s{simpl}.nc",
            country_shapes="resources/" + RDIR + "shapes/country_shapes.geojson",
            regions_onshore="resources/"
            + RDIR
            + "bus_regions/regions_onshore_elec_s{simpl}.geojson",
            regions_offshore="resources/"
            + RDIR
            + "bus_regions/regions_offshore_elec_s{simpl}.geojson",
            #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson",
            #using this line instead of the following will test updated gadm shapes for MA.
            #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
            #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
            gadm_shapes="resources/" + RDIR + "shapes/gadm_shapes.geojson",
            # busmap=ancient('resources/" + RDIR + "bus_regions/busmap_elec_s{simpl}.csv'),
            # custom_busmap=("data/custom_busmap_elec_s{simpl}_{clusters}.csv"
            #                if config["enable"].get("custom_busmap", False) else []),
            tech_costs=COSTS,
        output:
            network="networks/" + RDIR + "elec_s{simpl}_{clusters}_pre_augmentation.nc",
            regions_onshore="resources/"
            + RDIR
            + "bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
            regions_offshore="resources/"
            + RDIR
            + "bus_regions/regions_offshore_elec_s{simpl}_{clusters}.geojson",
            busmap="resources/"
            + RDIR
            + "bus_regions/busmap_elec_s{simpl}_{clusters}.csv",
            linemap="resources/"
            + RDIR
            + "bus_regions/linemap_elec_s{simpl}_{clusters}.csv",
        log:
            "logs/" + RDIR + "cluster_network/elec_s{simpl}_{clusters}.log",
        benchmark:
            "benchmarks/" + RDIR + "cluster_network/elec_s{simpl}_{clusters}"
        threads: 1
        resources:
            mem_mb=3000,
        script:
            "scripts/cluster_network.py"

    rule augmented_line_connections:
        input:
            tech_costs=COSTS,
            network="networks/" + RDIR + "elec_s{simpl}_{clusters}_pre_augmentation.nc",
            regions_onshore="resources/"
            + RDIR
            + "bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
            regions_offshore="resources/"
            + RDIR
            + "bus_regions/regions_offshore_elec_s{simpl}_{clusters}.geojson",
        output:
            network="networks/" + RDIR + "elec_s{simpl}_{clusters}.nc",
        log:
            "logs/" + RDIR + "augmented_line_connections/elec_s{simpl}_{clusters}.log",
        benchmark:
            "benchmarks/" + RDIR + "augmented_line_connections/elec_s{simpl}_{clusters}"
        threads: 1
        resources:
            mem_mb=3000,
        script:
            "scripts/augmented_line_connections.py"


if config["augmented_line_connection"].get("add_to_snakefile", False) == False:

    rule cluster_network:
        input:
            network="networks/" + RDIR + "elec_s{simpl}.nc",
            country_shapes="resources/" + RDIR + "shapes/country_shapes.geojson",
            regions_onshore="resources/"
            + RDIR
            + "bus_regions/regions_onshore_elec_s{simpl}.geojson",
            regions_offshore="resources/"
            + RDIR
            + "bus_regions/regions_offshore_elec_s{simpl}.geojson",
            #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson",
            #using this line instead of the following will test updated gadm shapes for MA.
            #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
            #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
            gadm_shapes="resources/" + RDIR + "shapes/gadm_shapes.geojson",
            # busmap=ancient('resources/" + RDIR + "bus_regions/busmap_elec_s{simpl}.csv'),
            # custom_busmap=("data/custom_busmap_elec_s{simpl}_{clusters}.csv"
            #                if config["enable"].get("custom_busmap", False) else []),
            tech_costs=COSTS,
        output:
            network="networks/" + RDIR + "elec_s{simpl}_{clusters}.nc",
            regions_onshore="resources/"
            + RDIR
            + "bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
            regions_offshore="resources/"
            + RDIR
            + "bus_regions/regions_offshore_elec_s{simpl}_{clusters}.geojson",
            busmap="resources/"
            + RDIR
            + "bus_regions/busmap_elec_s{simpl}_{clusters}.csv",
            linemap="resources/"
            + RDIR
            + "bus_regions/linemap_elec_s{simpl}_{clusters}.csv",
        log:
            "logs/" + RDIR + "cluster_network/elec_s{simpl}_{clusters}.log",
        benchmark:
            "benchmarks/" + RDIR + "cluster_network/elec_s{simpl}_{clusters}"
        threads: 1
        resources:
            mem_mb=3000,
        script:
            "scripts/cluster_network.py"


rule add_extra_components:
    input:
        network="networks/" + RDIR + "elec_s{simpl}_{clusters}.nc",
        tech_costs=COSTS,
    output:
        "networks/" + RDIR + "elec_s{simpl}_{clusters}_ec.nc",
    log:
        "logs/" + RDIR + "add_extra_components/elec_s{simpl}_{clusters}.log",
    benchmark:
        "benchmarks/" + RDIR + "add_extra_components/elec_s{simpl}_{clusters}_ec"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/add_extra_components.py"


rule prepare_network:
    input:
        "networks/" + RDIR + "elec_s{simpl}_{clusters}_ec.nc",
        tech_costs=COSTS,
    output:
        "networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    log:
        "logs/" + RDIR + "prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.log",
    benchmark:
        (
            "benchmarks/"
            + RDIR
            + "prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
        )
    threads: 1
    resources:
        mem_mb=4000,
    script:
        "scripts/prepare_network.py"


def memory(w):
    factor = 3.0
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)h$", o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)seg$", o, re.IGNORECASE)
        if m is not None:
            factor *= int(m.group(1)) / 8760
            break
    if w.clusters.endswith("m"):
        return int(factor * (18000 + 180 * int(w.clusters[:-1])))
    else:
        return int(factor * (10000 + 195 * int(w.clusters)))


if config["monte_carlo"]["options"].get("add_to_snakefile", False) == False:

    rule solve_network:
        input:
            "networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        output:
            "results/" + RDIR + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        log:
            solver=normpath(
                "logs/"
                + RDIR
                + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_solver.log"
            ),
            python="logs/"
            + RDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_python.log",
            memory="logs/"
            + RDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_memory.log",
        benchmark:
            (
                "benchmarks/"
                + RDIR
                + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
            )
        threads: 20
        resources:
            mem=memory,
        shadow:
            "shallow"
        script:
            "scripts/solve_network.py"


if config["monte_carlo"]["options"].get("add_to_snakefile", False) == True:

    rule monte_carlo:
        input:
            "networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        output:
            "networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}.nc",
        log:
            "logs/"
            + RDIR
            + "prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}.log",
        benchmark:
            (
                "benchmarks/"
                + RDIR
                + "prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}"
            )
        threads: 1
        resources:
            mem_mb=4000,
        script:
            "scripts/monte_carlo.py"

    rule solve_monte:
        input:
            expand(
                "networks/"
                + RDIR
                + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}.nc",
                **config["scenario"]
            ),

    rule solve_network:
        input:
            "networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}.nc",
        output:
            "results/"
            + RDIR
            + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}.nc",
        log:
            solver=normpath(
                "logs/"
                + RDIR
                + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}_solver.log"
            ),
            python="logs/"
            + RDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}_python.log",
            memory="logs/"
            + RDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}_memory.log",
        benchmark:
            (
                "benchmarks/"
                + RDIR
                + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}"
            )
        threads: 20
        resources:
            mem_mb=memory,
        shadow:
            "shallow"
        script:
            "scripts/solve_network.py"

    rule solve_all_networks_monte:
        input:
            expand(
                "results/"
                + RDIR
                + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}.nc",
                **config["scenario"]
            ),


def input_make_summary(w):
    # It's mildly hacky to include the separate costs input as first entry
    if w.ll.endswith("all"):
        ll = config["scenario"]["ll"]
        if len(w.ll) == 4:
            ll = [l for l in ll if l[0] == w.ll[0]]
    else:
        ll = w.ll
    return [COSTS] + expand(
        "results/" + RDIR + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        ll=ll,
        **{
            k: config["scenario"][k] if getattr(w, k) == "all" else getattr(w, k)
            for k in ["simpl", "clusters", "opts"]
        }
    )


rule make_summary:
    input:
        input_make_summary,
        tech_costs=COSTS,
    output:
        directory(
            "results/"
            + RDIR
            + "summaries/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}"
        ),
    log:
        "logs/"
        + RDIR
        + "make_summary/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}.log",
    script:
        "scripts/make_summary.py"


rule plot_summary:
    input:
        "results/"
        + RDIR
        + "summaries/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}",
    output:
        "results/"
        + RDIR
        + "plots/summary_{summary}_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}.{ext}",
    log:
        "logs/"
        + RDIR
        + "plot_summary/{summary}_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}_{ext}.log",
    script:
        "scripts/plot_summary.py"


rule plot_network:
    input:
        network="results/"
        + RDIR
        + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        africa_shape="resources/" + RDIR + "shapes/africa_shape.geojson",
        tech_costs=COSTS,
    output:
        only_map="results/"
        + RDIR
        + "plots/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}.{ext}",
        ext="results/"
        + RDIR
        + "plots/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}_ext.{ext}",
    log:
        "logs/"
        + RDIR
        + "plot_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}_{ext}.log",
    script:
        "scripts/plot_network.py"


rule build_test_configs:
    input:
        base_config="config.tutorial.yaml",
        update_file_list=[
            "test/config.tutorial_noprogress.yaml",
            "test/config.custom.yaml",
            "test/config.monte_carlo.yaml",
            "test/config.landlock.yaml",
        ],
    output:
        tmp_test_configs=[
            "test/tmp/config.tutorial_noprogress_tmp.yaml",
            "test/tmp/config.custom_tmp.yaml",
            "test/tmp/config.monte_carlo_tmp.yaml",
            "test/tmp/config.landlock_tmp.yaml",
        ],
    script:
        "scripts/build_test_configs.py"


rule make_statistics:
    output:
        stats="results/" + RDIR + "stats.csv",
    threads: 1
    script:
        "scripts/make_statistics.py"


rule run_scenario:
    input:
        default_config=DEFAULT_CONFIG,
        diff_config="configs/scenarios/config.{scenario_name}.yaml",
    output:
        touchfile=touch("results/{scenario_name}/scenario.done"),
        copyconfig="results/{scenario_name}/config.yaml",
    threads: 1
    resources:
        mem_mb=5000,
    run:
        from scripts.build_test_configs import create_test_config

        # Ensure the scenario name matches the name of the configuration
        create_test_config(
            input.diff_config,
            {"run": {"name": wildcards.scenario_name}},
            input.diff_config,
        )
        # merge the default config file with the difference
        create_test_config(input.default_config, input.diff_config, "config.yaml")
        os.system("snakemake -j all solve_all_networks --rerun-incomplete")
        os.system("snakemake -j1 make_statistics --force")
        copyfile("config.yaml", output.copyconfig)


rule run_all_scenarios:
    input:
        expand(
            "results/{scenario_name}/scenario.done",
            scenario_name=[
                c.stem.replace("config.", "")
                for c in Path("configs/scenarios").glob("config.*.yaml")
            ],
        ),
