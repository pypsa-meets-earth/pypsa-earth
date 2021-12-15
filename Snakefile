# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later

from os.path import normpath, exists, isdir
from shutil import copyfile

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

if not exists("config.yaml"):
    copyfile("config.default.yaml", "config.yaml")


configfile: "config.yaml"


COSTS = "data/costs.csv"
ATLITE_NPROCESSES = config["atlite"].get("nprocesses", 20)


wildcard_constraints:
    simpl="[a-zA-Z0-9]*|all",
    clusters="[0-9]+m?|all",
    ll="(v|c)([0-9\.]+|opt|all)|all",
    opts="[-+a-zA-Z0-9\.]*",


rule solve_all_networks:
    input: expand("results/networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc", **config['scenario'])


datafiles = [
        "resources/ssp2-2.6/2030/era5_2013/Africa.nc",
        "data/raw/copernicus/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        "data/raw/gebco/GEBCO_2021_TID.nc",
        "data/raw/eez/eez_v11.gpkg",
        # "data/raw/landcover",  # set as an explicit directory in the rule
        "data/raw/hydrobasins/hybas_lake_af_lev04_v1c.shp",
        "data/costs.csv",
]

if config.get('tutorial')==False:
    datafiles.extend(["cutouts/africa-2013-era5.nc"])
if config.get('tutorial')==True:
    datafiles.extend(["cutouts/africa-2013-era5-tutorial.nc"])

if config['enable'].get('retrieve_databundle', True):
    rule retrieve_databundle_light:
        output: #expand(directory('{file}') if isdir('{file}') else '{file}', file=datafiles)
            expand('{file}', file=datafiles),
            directory("data/raw/landcover")
        log: "logs/retrieve_databundle.log"
        script: 'scripts/retrieve_databundle_light.py'

if config['enable'].get('download_osm_data', True):
    rule download_osm_data:
        output:
            cables="data/raw/africa_all_raw_cables.geojson",
            generators="data/raw/africa_all_raw_generators.geojson",
            generators_csv="data/raw/africa_all_raw_generators.csv",
            lines="data/raw/africa_all_raw_lines.geojson",
            substations="data/raw/africa_all_raw_substations.geojson",
        log: "logs/download_osm_data.log"
        script: "scripts/osm_pbf_power_data_extractor.py"


rule clean_osm_data:
    input:
        cables="data/raw/africa_all_raw_cables.geojson",
        generators="data/raw/africa_all_raw_generators.geojson",
        lines="data/raw/africa_all_raw_lines.geojson",
        substations="data/raw/africa_all_raw_substations.geojson",
        country_shapes='resources/country_shapes.geojson',
        offshore_shapes='resources/offshore_shapes.geojson',
        africa_shape='resources/africa_shape.geojson',
    output:
        generators="data/clean/africa_all_generators.geojson",
        generators_csv="data/clean/africa_all_generators.csv",
        lines="data/clean/africa_all_lines.geojson",
        substations="data/clean/africa_all_substations.geojson",
    log: "logs/clean_osm_data.log"
    script: "scripts/osm_data_cleaning.py"


rule build_osm_network:
    input:
        generators="data/clean/africa_all_generators.geojson",
        lines="data/clean/africa_all_lines.geojson",
        substations="data/clean/africa_all_substations.geojson",
        country_shapes='resources/country_shapes.geojson',
    output:
        lines="data/base_network/africa_all_lines_build_network.csv",
        substations="data/base_network/africa_all_buses_build_network.csv",
    log: "logs/build_osm_network.log"
    script: "scripts/osm_built_network.py"


rule build_shapes:
    input:
        # naturalearth='data/bundle/naturalearth/ne_10m_admin_0_countries.shp',
        # eez='data/bundle/eez/World_EEZ_v8_2014.shp',
        # nuts3='data/bundle/NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp',
        # nuts3pop='data/bundle/nama_10r_3popgdp.tsv.gz',
        # nuts3gdp='data/bundle/nama_10r_3gdp.tsv.gz',
        eez='data/raw/eez/eez_v11.gpkg'
    output:
        country_shapes='resources/country_shapes.geojson',
        offshore_shapes='resources/offshore_shapes.geojson',
        africa_shape='resources/africa_shape.geojson',
        gadm_shapes='resources/gadm_shapes.geojson'
    log: "logs/build_shapes.log"
    threads: 1
    resources: mem=500
    script: "scripts/build_shapes.py"


rule base_network:
    input:
        osm_buses="data/base_network/africa_all_buses_build_network.csv",
        osm_lines="data/base_network/africa_all_lines_build_network.csv",
        country_shapes='resources/country_shapes.geojson',
        offshore_shapes='resources/offshore_shapes.geojson',
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
        country_shapes='resources/country_shapes.geojson',
        offshore_shapes='resources/offshore_shapes.geojson',
        base_network="networks/base.nc",
        gadm_shapes="resources/gadm_shapes.geojson"
    output:
        regions_onshore="resources/regions_onshore.geojson",
        regions_offshore="resources/regions_offshore.geojson"
    log: "logs/build_bus_regions.log"
    threads: 1
    resources: mem=1000
    script: "scripts/build_bus_regions.py"


if config['enable'].get('build_cutout', False):
    rule build_cutout:
        input:
            regions_onshore="resources/regions_onshore.geojson",
            regions_offshore="resources/regions_offshore.geojson"
        output: "cutouts/{cutout}.nc"
        log: "logs/build_cutout/{cutout}.log"
        benchmark: "benchmarks/build_cutout_{cutout}"
        threads: ATLITE_NPROCESSES
        resources: mem=ATLITE_NPROCESSES * 1000
        script: "scripts/build_cutout.py"


if config['enable'].get('build_natura_raster', False):
    rule build_natura_raster:
        input:
            shapefiles_land="data/raw/landcover",
            cutouts=expand("cutouts/{cutouts}.nc", **config['atlite'])
        output: "resources/natura.tiff"
        log: "logs/build_natura_raster.log"
        script: "scripts/build_natura_raster.py"


rule build_renewable_profiles:
    input:
        base_network="networks/base.nc",
        natura="resources/natura.tiff",
        copernicus="data/raw/copernicus/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        gebco="data/raw/gebco/GEBCO_2021_TID.nc",
        country_shapes='resources/country_shapes.geojson',
        offshore_shapes='resources/offshore_shapes.geojson',
        regions=lambda w: ("resources/regions_onshore.geojson"
                                   if w.technology in ('onwind', 'solar', "hydro")
                                   else "resources/regions_offshore.geojson"),
        cutout=lambda w: "cutouts/" + config["renewable"][w.technology]['cutout'] + ".nc",
    output: profile="resources/profile_{technology}.nc",
    log: "logs/build_renewable_profile_{technology}.log"
    benchmark: "benchmarks/build_renewable_profiles_{technology}"
    threads: ATLITE_NPROCESSES
    resources: mem=ATLITE_NPROCESSES * 5000
    script: "scripts/build_renewable_profiles.py"


rule build_powerplants:
    input:
        base_network="networks/base.nc",
        pm_config="configs/powerplantmatching_config.yaml",
        custom_powerplants="data/clean/africa_all_generators.csv"
    output: "resources/powerplants.csv"
    log: "logs/build_powerplants.log"
    threads: 1
    resources: mem=500
    script: "scripts/build_powerplants.py"


rule add_electricity:
    input:
        base_network='networks/base.nc',
        tech_costs=COSTS,
        regions="resources/regions_onshore.geojson",
        powerplants='resources/powerplants.csv',
        # hydro_capacities='data/bundle/hydro_capacities.csv',
        # geth_hydro_capacities='data/geth2015_hydro_capacities.csv',
        load='resources/ssp2-2.6/2050/era5_2013/Africa.nc',
        gadm_shapes='resources/gadm_shapes.geojson',
        **{f"profile_{tech}": f"resources/profile_{tech}.nc"
            for tech in config['renewable']}
    output: "networks/elec.nc"
    log: "logs/add_electricity.log"
    benchmark: "benchmarks/add_electricity"
    threads: 1
    resources: mem=3000
    script: "scripts/add_electricity.py"


rule simplify_network:
    input:
        network='networks/elec.nc',
        tech_costs=COSTS,
        regions_onshore="resources/regions_onshore.geojson",
        regions_offshore="resources/regions_offshore.geojson"
    output:
        network='networks/elec_s{simpl}.nc',
        regions_onshore="resources/regions_onshore_elec_s{simpl}.geojson",
        regions_offshore="resources/regions_offshore_elec_s{simpl}.geojson",
        busmap='resources/busmap_elec_s{simpl}.csv'
    log: "logs/simplify_network/elec_s{simpl}.log"
    benchmark: "benchmarks/simplify_network/elec_s{simpl}"
    threads: 1
    resources: mem=4000
    script: "scripts/simplify_network.py"


rule cluster_network:
    input:
        network='networks/elec_s{simpl}.nc',
        raw_onshore_busregion="resources/regions_onshore.geojson",
        regions_onshore="resources/regions_onshore_elec_s{simpl}.geojson",
        regions_offshore="resources/regions_offshore_elec_s{simpl}.geojson",
        # busmap=ancient('resources/busmap_elec_s{simpl}.csv'),
        # custom_busmap=("data/custom_busmap_elec_s{simpl}_{clusters}.csv"
        #                if config["enable"].get("custom_busmap", False) else []),
        tech_costs=COSTS
    output:
        network='networks/elec_s{simpl}_{clusters}.nc',
        regions_onshore="resources/regions_onshore_elec_s{simpl}_{clusters}.geojson",
        regions_offshore="resources/regions_offshore_elec_s{simpl}_{clusters}.geojson",
        busmap="resources/busmap_elec_s{simpl}_{clusters}.csv",
        linemap="resources/linemap_elec_s{simpl}_{clusters}.csv"
    log: "logs/cluster_network/elec_s{simpl}_{clusters}.log"
    benchmark: "benchmarks/cluster_network/elec_s{simpl}_{clusters}"
    threads: 1
    resources: mem=3000
    script: "scripts/cluster_network.py"


rule augmented_line_connections:
    input:
        tech_costs=COSTS,
        network='networks/elec_s{simpl}_{clusters}.nc',
        regions_onshore="resources/regions_onshore_elec_s{simpl}_{clusters}.geojson",
        regions_offshore="resources/regions_offshore_elec_s{simpl}_{clusters}.geojson",
    output:
        network='networks/elec_s{simpl}_{clusters}.nc',
        network_unmeshed='networks/elec_s{simpl}_{clusters}_no_augmented_connections.nc',
        regions_onshore="resources/regions_onshore_elec_s{simpl}_{clusters}.geojson",
        regions_offshore="resources/regions_offshore_elec_s{simpl}_{clusters}.geojson",
    log: "logs/augmented_line_connections/elec_s{simpl}_{clusters}.log"
    benchmark: "benchmarks/augmented_line_connections/elec_s{simpl}_{clusters}"
    threads: 1
    resources: mem=3000
    script: "scripts/augmented_line_connections.py"


rule add_extra_components:
    input:
        network='networks/elec_s{simpl}_{clusters}.nc',
        tech_costs=COSTS,
    output: 'networks/elec_s{simpl}_{clusters}_ec.nc'
    log: "logs/add_extra_components/elec_s{simpl}_{clusters}.log"
    benchmark: "benchmarks/add_extra_components/elec_s{simpl}_{clusters}_ec"
    threads: 1
    resources: mem=3000
    script: "scripts/add_extra_components.py"


rule prepare_network:
    input: 'networks/elec_s{simpl}_{clusters}_ec.nc', tech_costs=COSTS
    output: 'networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc'
    log: "logs/prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.log"
    benchmark: "benchmarks/prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
    threads: 1
    resources: mem=4000
    script: "scripts/prepare_network.py"


def memory(w):
    factor = 3.
    for o in w.opts.split('-'):
        m = re.match(r'^(\d+)h$', o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    for o in w.opts.split('-'):
        m = re.match(r'^(\d+)seg$', o, re.IGNORECASE)
        if m is not None:
            factor *= int(m.group(1)) / 8760
            break
    if w.clusters.endswith('m'):
        return int(factor * (18000 + 180 * int(w.clusters[:-1])))
    else:
        return int(factor * (10000 + 195 * int(w.clusters)))


rule solve_network:
    input: "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"
    output: "results/networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"
    log:
        solver=normpath("logs/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_solver.log"),
        python="logs/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_python.log",
        memory="logs/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_memory.log"
    benchmark: "benchmarks/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
    threads: 20
    resources: mem=memory
    shadow: "shallow"
    script: "scripts/solve_network.py"
