# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later

from os.path import normpath, exists
from shutil import copyfile

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

if not exists("config.yaml"):
    copyfile("config.default.yaml", "config.yaml")


configfile: "config.yaml"


COSTS = "data/costs.csv"
ATLITE_NPROCESSES = config["atlite"].get("nprocesses", 4)


wildcard_constraints:
    simpl="[a-zA-Z0-9]*|all",
    clusters="[0-9]+m?|all",
    ll="(v|c)([0-9\.]+|opt|all)|all",
    opts="[-+a-zA-Z0-9\.]*",


if config['enable'].get('build_cutout', False):
    rule build_cutout:
        input:
            regions_onshore="data/shapes/regions_onshore.geojson",
            regions_offshore="data/shapes/regions_offshore.geojson"
        output: "cutouts/{cutout}.nc"
        log: "logs/build_cutout/{cutout}.log"
        benchmark: "benchmarks/build_cutout_{cutout}"
        threads: ATLITE_NPROCESSES
        resources: mem=ATLITE_NPROCESSES * 1000
        script: "scripts/build_cutout.py"


# if config['enable'].get('build_natura_raster', False):
#     rule build_natura_raster:
#         input:
#             natura="data/bundle/natura/Natura2000_end2015.shp",
#             cutouts=expand("cutouts/{cutouts}.nc", **config['atlite'])
#         output: "resources/natura.tiff"
#         log: "logs/build_natura_raster.log"
#         script: "scripts/build_natura_raster.py"


rule base_network:
    input:
        osm_buses="data/base_network/africa_all_buses_build_network.csv",
        osm_lines="data/base_network/africa_all_lines_build_network.csv",
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
        # country_shapes='resources/country_shapes.geojson',
        # offshore_shapes='resources/offshore_shapes.geojson',
        # europe_shape='resources/europe_shape.geojson'
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


rule add_electricity:
    input:
        base_network='networks/base.nc',
        tech_costs=COSTS,
        # regions="resources/regions_onshore.geojson",
        # powerplants='resources/powerplants.csv',
        # hydro_capacities='data/bundle/hydro_capacities.csv',
        # geth_hydro_capacities='data/geth2015_hydro_capacities.csv',
        # load='resources/load.csv',
        # nuts3_shapes='resources/nuts3_shapes.geojson',
        # **{f"profile_{tech}": f"resources/profile_{tech}.nc"
        #    for tech in config['renewable']}
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
        # regions_onshore="resources/regions_onshore.geojson",
        # regions_offshore="resources/regions_offshore.geojson"
    output:
        network='networks/elec_s{simpl}.nc',
        # regions_onshore="resources/regions_onshore_elec_s{simpl}.geojson",
        # regions_offshore="resources/regions_offshore_elec_s{simpl}.geojson",
        # busmap='resources/busmap_elec_s{simpl}.csv'
    log: "logs/simplify_network/elec_s{simpl}.log"
    benchmark: "benchmarks/simplify_network/elec_s{simpl}"
    threads: 1
    resources: mem=4000
    script: "scripts/simplify_network.py"


rule cluster_network:
    input:
        # network='networks/elec_s.nc',
        # regions_onshore="resources/regions_onshore_elec_s.geojson",
        # regions_offshore="resources/regions_offshore_elec_s.geojson",
        # busmap=ancient('resources/busmap_elec_s.csv'),
        # custom_busmap=("data/custom_busmap_elec_s_1000.csv"
        #                if config["enable"].get("custom_busmap", False) else []),
        # tech_costs=COSTS
    output:
        # network='networks/elec_s_1000.nc',
        # regions_onshore="resources/regions_onshore_elec_s_1000.geojson",
        # regions_offshore="resources/regions_offshore_elec_s_1000.geojson",
        # busmap="resources/busmap_elec_s_1000.csv",
        # linemap="resources/linemap_elec_s_1000.csv"
    log: "logs/cluster_network/elec_s_1000.log"
    benchmark: "benchmarks/cluster_network/elec_s_1000"
    threads: 1
    resources: mem=3000
    script: "scripts/cluster_network.py"

# rule build_shapes:
#     input:
#         naturalearth='data/bundle/naturalearth/ne_10m_admin_0_countries.shp',
#         eez='data/bundle/eez/World_EEZ_v8_2014.shp',
#         nuts3='data/bundle/NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp',
#         nuts3pop='data/bundle/nama_10r_3popgdp.tsv.gz',
#         nuts3gdp='data/bundle/nama_10r_3gdp.tsv.gz',
#     output:
#         country_shapes='resources/country_shapes.geojson',
#         offshore_shapes='resources/offshore_shapes.geojson',
#         europe_shape='resources/europe_shape.geojson',
#         nuts3_shapes='resources/nuts3_shapes.geojson'
#     log: "logs/build_shapes.log"
#     threads: 1
#     resources: mem=500
#     script: "scripts/build_shapes.py"



