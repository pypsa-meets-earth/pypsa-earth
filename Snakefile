# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later

import sys

sys.path.append("./scripts")

from os.path import normpath, exists, isdir
from shutil import copyfile

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

from scripts.download_osm_data import create_country_list
from scripts.add_electricity import get_load_paths_gegis

HTTP = HTTPRemoteProvider()

if "config" not in globals() or not config:  # skip when used as sub-workflow
    if not exists("config.yaml"):
        copyfile("config.default.yaml", "config.yaml")

    configfile: "config.yaml"


configfile: "configs/bundle_config.yaml"


# convert country list according to the desired region
config["countries"] = create_country_list(config["countries"])

load_data_paths = get_load_paths_gegis("data", config)
COSTS = "data/costs.csv"
ATLITE_NPROCESSES = config["atlite"].get("nprocesses", 20)


wildcard_constraints:
    simpl="[a-zA-Z0-9]*|all",
    clusters="[0-9]+m?|all",
    ll="(v|c)([0-9\.]+|opt|all)|all",
    opts="[-+a-zA-Z0-9\.]*",


rule run_test:
    run:
        shell("cp test/config.test1.yaml config.yaml")
        shell("snakemake --cores all solve_all_networks --forceall")


rule solve_all_networks:
    input:
        expand(
            "results/networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            **config["scenario"]
        ),


rule plot_all_p_nom:
    input:
        expand(
            "results/plots/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_p_nom.{ext}",
            **config["scenario"],
            ext=["png", "pdf"]
        ),


rule make_all_summaries:
    input:
        expand(
            "results/summaries/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}",
            **config["scenario"],
            country=["all"] + config["countries"]
        ),


rule plot_all_summaries:
    input:
        expand(
            "results/plots/summary_{summary}_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}.{ext}",
            summary=["energy", "costs"],
            **config["scenario"],
            country=["all"] + config["countries"],
            ext=["png", "pdf"]
        ),


def datafiles_retrivedatabundle(config):
    listoutputs = [
        dvalue["output"]
        for (dname, dvalue) in config["databundles"].items()
        if config.get("tutorial", False) == dvalue.get("tutorial", False)
    ]
    unique_outputs = set(
        (
            [
                inneroutput
                for output in listoutputs
                for inneroutput in output
                if "*" not in inneroutput
                or inneroutput.endswith("/")  # exclude directories
            ]
        )
    )

    # when option build_natura_raster is enabled, remove natura.tiff from the outputs
    if config["enable"].get("build_natura_raster", False):
        unique_outputs = [
            output for output in unique_outputs if "natura.tiff" not in output
        ]

    # when option build_cutout is enabled, remove cutouts from the outputs
    if config["enable"].get("build_cutout", False):
        unique_outputs = [
            output for output in unique_outputs if "cutouts/" not in output
        ]

    return unique_outputs


if config["enable"].get("retrieve_databundle", True):

    rule retrieve_databundle_light:
        output:  #expand(directory('{file}') if isdir('{file}') else '{file}', file=datafiles)
            expand("{file}", file=datafiles_retrivedatabundle(config)),
            directory("data/landcover"),
        log:
            "logs/retrieve_databundle.log",
        script:
            "scripts/retrieve_databundle_light.py"


if config["enable"].get("download_osm_data", True):

    rule download_osm_data:
        output:
            cables="resources/osm/raw/africa_all_raw_cables.geojson",
            generators="resources/osm/raw/africa_all_raw_generators.geojson",
            generators_csv="resources/osm/raw/africa_all_raw_generators.csv",
            lines="resources/osm/raw/africa_all_raw_lines.geojson",
            substations="resources/osm/raw/africa_all_raw_substations.geojson",
        log:
            "logs/download_osm_data.log",
        script:
            "scripts/download_osm_data.py"


rule clean_osm_data:
    input:
        cables="resources/osm/raw/africa_all_raw_cables.geojson",
        generators="resources/osm/raw/africa_all_raw_generators.geojson",
        lines="resources/osm/raw/africa_all_raw_lines.geojson",
        substations="resources/osm/raw/africa_all_raw_substations.geojson",
        country_shapes="resources/shapes/country_shapes.geojson",
        offshore_shapes="resources/shapes/offshore_shapes.geojson",
        africa_shape="resources/shapes/africa_shape.geojson",
    output:
        generators="resources/osm/clean/africa_all_generators.geojson",
        generators_csv="resources/osm/clean/africa_all_generators.csv",
        lines="resources/osm/clean/africa_all_lines.geojson",
        substations="resources/osm/clean/africa_all_substations.geojson",
    log:
        "logs/clean_osm_data.log",
    script:
        "scripts/clean_osm_data.py"


rule build_osm_network:
    input:
        generators="resources/osm/clean/africa_all_generators.geojson",
        lines="resources/osm/clean/africa_all_lines.geojson",
        substations="resources/osm/clean/africa_all_substations.geojson",
        country_shapes="resources/shapes/country_shapes.geojson",
    output:
        lines="resources/base_network/africa_all_lines_build_network.csv",
        converters="resources/base_network/africa_all_converters_build_network.csv",
        transformers="resources/base_network/africa_all_transformers_build_network.csv",
        substations="resources/base_network/africa_all_buses_build_network.csv",
    log:
        "logs/build_osm_network.log",
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
        country_shapes="resources/shapes/country_shapes.geojson",
        offshore_shapes="resources/shapes/offshore_shapes.geojson",
        africa_shape="resources/shapes/africa_shape.geojson",
        gadm_shapes="resources/shapes/gadm_shapes.geojson",
    log:
        "logs/build_shapes.log",
    threads: 1
    resources:
        mem=500,
    script:
        "scripts/build_shapes.py"


rule base_network:
    input:
        osm_buses="resources/base_network/africa_all_buses_build_network.csv",
        osm_lines="resources/base_network/africa_all_lines_build_network.csv",
        osm_converters="resources/base_network/africa_all_converters_build_network.csv",
        osm_transformers="resources/base_network/africa_all_transformers_build_network.csv",
        country_shapes="resources/shapes/country_shapes.geojson",
        offshore_shapes="resources/shapes/offshore_shapes.geojson",
        # osm_buses='data/osm/africa_all_buses_clean.csv',
        # osm_lines='data/osm/africa_all_lines_clean.csv',
        # eg_buses='data/entsoegridkit/buses.csv',
        # eg_lines='data/entsoegridkit/lines.csv',
        # eg_links='data/entsoegridkit/links.csv',
        # eg_converters='data/entsoegridkit/converters.csv',
        # eg_transformers='data/entsoegridkit/transformers.csv',
        # parameter_corrections='data/parameter_corrections.yaml',
        # links_p_nom='data/links_p_nom.csv',
        # links_tyndp='data/links_tyndp.csv',
        # africa_shape='resources/africa_shape.geojson'
    output:
        "networks/base.nc",
    log:
        "logs/base_network.log",
    benchmark:
        "benchmarks/base_network"
    threads: 1
    resources:
        mem=500,
    script:
        "scripts/base_network.py"


rule build_bus_regions:
    input:
        country_shapes="resources/shapes/country_shapes.geojson",
        offshore_shapes="resources/shapes/offshore_shapes.geojson",
        base_network="networks/base.nc",
        #gadm_shapes="resources/shapes/MAR2.geojson",
        #using this line instead of the following will test updated gadm shapes for MA.
        #To use: downlaod file from the google drive and place it in resources/shapes/
        #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
        gadm_shapes="resources/shapes/gadm_shapes.geojson",
    output:
        regions_onshore="resources/bus_regions/regions_onshore.geojson",
        regions_offshore="resources/bus_regions/regions_offshore.geojson",
    log:
        "logs/build_bus_regions.log",
    threads: 1
    resources:
        mem=1000,
    script:
        "scripts/build_bus_regions.py"


if config["enable"].get("build_cutout", False):

    rule build_cutout:
        input:
            onshore_shapes="resources/shapes/country_shapes.geojson",
            offshore_shapes="resources/shapes/offshore_shapes.geojson",
        output:
            "cutouts/{cutout}.nc",
        log:
            "logs/build_cutout/{cutout}.log",
        benchmark:
            "benchmarks/build_cutout_{cutout}"
        threads: ATLITE_NPROCESSES
        resources:
            mem=ATLITE_NPROCESSES * 1000,
        script:
            "scripts/build_cutout.py"


if config["enable"].get("build_natura_raster", False):

    rule build_natura_raster:
        input:
            shapefiles_land="data/landcover",
            cutouts=expand("cutouts/{cutouts}.nc", **config["atlite"]),
        output:
            "resources/natura.tiff",
        log:
            "logs/build_natura_raster.log",
        script:
            "scripts/build_natura_raster.py"


if not config["enable"].get("build_natura_raster", False):

    rule copy_defaultnatura_tiff:
        input:
            "data/natura.tiff",
        output:
            "resources/natura.tiff",
        run:
            import shutil

            shutil.copyfile(input[0], output[0])


rule build_renewable_profiles:
    input:
        base_network="networks/base.nc",
        natura="resources/natura.tiff",
        copernicus="data/copernicus/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        gebco="data/gebco/GEBCO_2021_TID.nc",
        country_shapes="resources/shapes/country_shapes.geojson",
        offshore_shapes="resources/shapes/offshore_shapes.geojson",
        hydro_capacities="data/hydro_capacities.csv",
        eia_hydro_generation="data/eia_hydro_annual_generation.csv",
        powerplants="resources/powerplants.csv",
        regions=lambda w: (
            "resources/bus_regions/regions_onshore.geojson"
            if w.technology in ("onwind", "solar", "hydro")
            else "resources/bus_regions/regions_offshore.geojson"
        ),
        cutout=lambda w: "cutouts/"
        + config["renewable"][w.technology]["cutout"]
        + ".nc",
    output:
        profile="resources/renewable_profiles/profile_{technology}.nc",
    log:
        "logs/build_renewable_profile_{technology}.log",
    benchmark:
        "benchmarks/build_renewable_profiles_{technology}"
    threads: ATLITE_NPROCESSES
    resources:
        mem=ATLITE_NPROCESSES * 5000,
    script:
        "scripts/build_renewable_profiles.py"


rule build_powerplants:
    input:
        base_network="networks/base.nc",
        pm_config="configs/powerplantmatching_config.yaml",
        custom_powerplants="data/custom_powerplants.csv",
        osm_powerplants="resources/osm/clean/africa_all_generators.csv",
        #gadm_shapes="resources/shapes/MAR2.geojson",
        #using this line instead of the following will test updated gadm shapes for MA.
        #To use: downlaod file from the google drive and place it in resources/shapes/
        #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
        gadm_shapes="resources/shapes/gadm_shapes.geojson",
    output:
        powerplants="resources/powerplants.csv",
        powerplants_osm2pm="resources/powerplants_osm2pm.csv",
    log:
        "logs/build_powerplants.log",
    threads: 1
    resources:
        mem=500,
    script:
        "scripts/build_powerplants.py"


rule add_electricity:
    input:
        **{
            f"profile_{tech}": f"resources/renewable_profiles/profile_{tech}.nc"
            for tech in config["renewable"]
        },
        **{
            f"conventional_{carrier}_{attr}": fn
            for carrier, d in config.get("conventional", {None: {}}).items()
            for attr, fn in d.items()
            if str(fn).startswith("data/")
        },
        base_network="networks/base.nc",
        tech_costs=COSTS,
        regions="resources/bus_regions/regions_onshore.geojson",
        powerplants="resources/powerplants.csv",
        load=load_data_paths,
        #gadm_shapes="resources/shapes/MAR2.geojson", 
        #using this line instead of the following will test updated gadm shapes for MA.
        #To use: downlaod file from the google drive and place it in resources/shapes/
        #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
        gadm_shapes="resources/shapes/gadm_shapes.geojson",
        hydro_capacities="data/hydro_capacities.csv",
    output:
        "networks/elec.nc",
    log:
        "logs/add_electricity.log",
    benchmark:
        "benchmarks/add_electricity"
    threads: 1
    resources:
        mem=3000,
    script:
        "scripts/add_electricity.py"


rule simplify_network:
    input:
        network="networks/elec.nc",
        tech_costs=COSTS,
        regions_onshore="resources/bus_regions/regions_onshore.geojson",
        regions_offshore="resources/bus_regions/regions_offshore.geojson",
    output:
        network="networks/elec_s{simpl}.nc",
        regions_onshore="resources/bus_regions/regions_onshore_elec_s{simpl}.geojson",
        regions_offshore="resources/bus_regions/regions_offshore_elec_s{simpl}.geojson",
        busmap="resources/busmap_elec_s{simpl}.csv",
        connection_costs="resources/connection_costs_s{simpl}.csv",
    log:
        "logs/simplify_network/elec_s{simpl}.log",
    benchmark:
        "benchmarks/simplify_network/elec_s{simpl}"
    threads: 1
    resources:
        mem=4000,
    script:
        "scripts/simplify_network.py"


if config["augmented_line_connection"].get("add_to_snakefile", False) == True:

    rule cluster_network:
        input:
            network="networks/elec_s{simpl}.nc",
            country_shapes="resources/shapes/country_shapes.geojson",
            regions_onshore="resources/bus_regions/regions_onshore_elec_s{simpl}.geojson",
            regions_offshore="resources/bus_regions/regions_offshore_elec_s{simpl}.geojson",
            #gadm_shapes="resources/shapes/MAR2.geojson",
            #using this line instead of the following will test updated gadm shapes for MA.
            #To use: downlaod file from the google drive and place it in resources/shapes/
            #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
            gadm_shapes="resources/shapes/gadm_shapes.geojson",
            # busmap=ancient('resources/busmap_elec_s{simpl}.csv'),
            # custom_busmap=("data/custom_busmap_elec_s{simpl}_{clusters}.csv"
            #                if config["enable"].get("custom_busmap", False) else []),
            tech_costs=COSTS,
        output:
            network="networks/elec_s{simpl}_{clusters}_pre_augmentation.nc",
            regions_onshore="resources/bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
            regions_offshore="resources/bus_regions/regions_offshore_elec_s{simpl}_{clusters}.geojson",
            busmap="resources/busmap_elec_s{simpl}_{clusters}.csv",
            linemap="resources/linemap_elec_s{simpl}_{clusters}.csv",
        log:
            "logs/cluster_network/elec_s{simpl}_{clusters}.log",
        benchmark:
            "benchmarks/cluster_network/elec_s{simpl}_{clusters}"
        threads: 1
        resources:
            mem=3000,
        script:
            "scripts/cluster_network.py"

    rule augmented_line_connections:
        input:
            tech_costs=COSTS,
            network="networks/elec_s{simpl}_{clusters}_pre_augmentation.nc",
            regions_onshore="resources/bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
            regions_offshore="resources/bus_regions/regions_offshore_elec_s{simpl}_{clusters}.geojson",
        output:
            network="networks/elec_s{simpl}_{clusters}.nc",
        log:
            "logs/augmented_line_connections/elec_s{simpl}_{clusters}.log",
        benchmark:
            "benchmarks/augmented_line_connections/elec_s{simpl}_{clusters}"
        threads: 1
        resources:
            mem=3000,
        script:
            "scripts/augmented_line_connections.py"


if config["augmented_line_connection"].get("add_to_snakefile", False) == False:

    rule cluster_network:
        input:
            network="networks/elec_s{simpl}.nc",
            country_shapes="resources/shapes/country_shapes.geojson",
            regions_onshore="resources/bus_regions/regions_onshore_elec_s{simpl}.geojson",
            regions_offshore="resources/bus_regions/regions_offshore_elec_s{simpl}.geojson",
            # busmap=ancient('resources/busmap_elec_s{simpl}.csv'),
            # custom_busmap=("data/custom_busmap_elec_s{simpl}_{clusters}.csv"
            #                if config["enable"].get("custom_busmap", False) else []),
            tech_costs=COSTS,
        output:
            network="networks/elec_s{simpl}_{clusters}.nc",
            regions_onshore="resources/bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
            regions_offshore="resources/bus_regions/regions_offshore_elec_s{simpl}_{clusters}.geojson",
            busmap="resources/busmap_elec_s{simpl}_{clusters}.csv",
            linemap="resources/linemap_elec_s{simpl}_{clusters}.csv",
        log:
            "logs/cluster_network/elec_s{simpl}_{clusters}.log",
        benchmark:
            "benchmarks/cluster_network/elec_s{simpl}_{clusters}"
        threads: 1
        resources:
            mem=3000,
        script:
            "scripts/cluster_network.py"


rule add_extra_components:
    input:
        network="networks/elec_s{simpl}_{clusters}.nc",
        tech_costs=COSTS,
    output:
        "networks/elec_s{simpl}_{clusters}_ec.nc",
    log:
        "logs/add_extra_components/elec_s{simpl}_{clusters}.log",
    benchmark:
        "benchmarks/add_extra_components/elec_s{simpl}_{clusters}_ec"
    threads: 1
    resources:
        mem=3000,
    script:
        "scripts/add_extra_components.py"


rule prepare_network:
    input:
        "networks/elec_s{simpl}_{clusters}_ec.nc",
        tech_costs=COSTS,
    output:
        "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    log:
        "logs/prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.log",
    benchmark:
        "benchmarks/prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
    threads: 1
    resources:
        mem=4000,
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


rule solve_network:
    input:
        "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    output:
        "results/networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    log:
        solver=normpath(
            "logs/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_solver.log"
        ),
        python="logs/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_python.log",
        memory="logs/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_memory.log",
    benchmark:
        "benchmarks/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
    threads: 20
    resources:
        mem=memory,
    shadow:
        "shallow"
    script:
        "scripts/solve_network.py"


def input_make_summary(w):
    # It's mildly hacky to include the separate costs input as first entry
    if w.ll.endswith("all"):
        ll = config["scenario"]["ll"]
        if len(w.ll) == 4:
            ll = [l for l in ll if l[0] == w.ll[0]]
    else:
        ll = w.ll
    return [COSTS] + expand(
        "results/networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
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
            "results/summaries/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}"
        ),
    log:
        "logs/make_summary/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}.log",
    script:
        "scripts/make_summary.py"


rule plot_summary:
    input:
        "results/summaries/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}",
    output:
        "results/plots/summary_{summary}_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}.{ext}",
    log:
        "logs/plot_summary/{summary}_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}_{ext}.log",
    script:
        "scripts/plot_summary.py"


rule plot_network:
    input:
        network="results/networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        africa_shape="resources/shapes/africa_shape.geojson",
        tech_costs=COSTS,
    output:
        only_map="results/plots/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}.{ext}",
        ext="results/plots/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}_ext.{ext}",
    log:
        "logs/plot_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}_{ext}.log",
    script:
        "scripts/plot_network.py"
