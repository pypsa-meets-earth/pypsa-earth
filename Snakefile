
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

configfile: "config.yaml"

SDIR = config['summary_dir'] + '/' + config['run']
RDIR = config['results_dir'] + config['run']
CDIR = config['costs_dir']


rule prepare_sector_network:
    input:
        network='networks/elec_s{simpl}_{clusters}.nc',
        costs=CDIR + "costs_2030.csv",
        h2_cavern="data/hydrogen_salt_cavern_potentials.csv",
        energy_totals_name="resources/energy_totals.csv",
        traffic_data_KFZ = "data/emobility/KFZ__count",
        traffic_data_Pkw = "data/emobility/Pkw__count",

        transport_name='resources/transport_data.csv',

    output: RDIR + '/prenetworks/elec_s{simpl}_{clusters}.nc'

    script: "scripts/prepare_sector_network.py"
<<<<<<< HEAD
<<<<<<< HEAD

rule prepare_transport_data:
    input:
        traffic_data_KFZ = "data/emobility/KFZ__count",
        network='networks/elec_s{simpl}_{clusters}.nc',
        # Get pop layouts from Morocco (dummy)
        clustered_pop_layout_dummy="resources/pop_layout_elec_s{simpl}_dummy.csv",

    output: 

    script: "scripts/prepare_transport_data.py"

rule calculate_dummy_pop_layout:
    input:
        network='networks/elec_s{simpl}_{clusters}.nc',

        # Get pop layouts from Europe (update to Morocco/Africa layout)
        clustered_pop_layout="resources/pop_layout_elec_s{simpl}_37.csv",
        #simplified_pop_layout="resources/pop_layout_elec_s{simpl}.csv",

    output: clustered_pop_layout_dummy="resources/pop_layout_elec_s{simpl}_dummy.csv",

    script: "scripts/calculate_dummy_pop_layout.py" 

rule build_population_layouts:
    input:
        nuts3_shapes='resources/gadm_shapes.geojson',
        urban_percent="data/urban_percent.csv"
    output:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc"
    resources: mem_mb=20000
    benchmark: "benchmarks/build_population_layouts"
    threads: 8
    script: "scripts/build_population_layouts.py"