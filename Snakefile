
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
        nodal_energy_totals='resources/nodal_energy_totals.csv',
        transport='resources/transport.csv',
        avail_profile='resources/avail_profile.csv',
        dsm_profile='resources/dsm_profile.csv',
        nodal_transport_data='resources/nodal_transport_data.csv',
        

    output: RDIR + '/prenetworks/elec_s{simpl}_{clusters}.nc'

    script: "scripts/prepare_sector_network.py"

rule prepare_transport_data:
    input:
        network='networks/elec_s{simpl}_{clusters}.nc',
        energy_totals_name='resources/energy_totals.csv',
        traffic_data_KFZ = "data/emobility/KFZ__count",
        traffic_data_Pkw = "data/emobility/Pkw__count",
        transport_name='resources/transport_data.csv',
        # Get pop layouts from Morocco (dummy)
        clustered_pop_layout_dummy="resources/pop_layout_elec_s{simpl}_dummy.csv",
        temp_air_total="resources/temp_air_total_elec_s{simpl}_37.nc",

    output: 
        nodal_energy_totals='resources/nodal_energy_totals.csv',
        transport='resources/transport.csv',
        avail_profile='resources/avail_profile.csv',
        dsm_profile='resources/dsm_profile.csv',
        nodal_transport_data='resources/nodal_transport_data.csv',

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
    
    
rule build_temperature_profiles:
        input:
            pop_layout_total="resources/pop_layout_total.nc",
            pop_layout_urban="resources/pop_layout_urban.nc",
            pop_layout_rural="resources/pop_layout_rural.nc",
            regions_onshore="resources/regions_onshore_elec_s{simpl}_{clusters}.geojson"
        output:
            temp_soil_total="resources/temp_soil_total_elec_s{simpl}_{clusters}.nc",
            temp_soil_rural="resources/temp_soil_rural_elec_s{simpl}_{clusters}.nc",
            temp_soil_urban="resources/temp_soil_urban_elec_s{simpl}_{clusters}.nc",
            temp_air_total="resources/temp_air_total_elec_s{simpl}_{clusters}.nc",
            temp_air_rural="resources/temp_air_rural_elec_s{simpl}_{clusters}.nc",
            temp_air_urban="resources/temp_air_urban_elec_s{simpl}_{clusters}.nc"
        resources: mem_mb=20000
        benchmark: "benchmarks/build_temperature_profiles/s{simpl}_{clusters}"
        script: "scripts/build_temperature_profiles.py"

    
rule solve_network:
        input:
            overrides="data/override_component_attrs",
            network=RDIR + "/prenetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc",
            costs=CDIR + "costs_{planning_horizons}.csv",
            config=SDIR + '/configs/config.yaml'
        output: RDIR + "/postnetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc"
        shadow: "shallow"
        log:
            solver=RDIR + "/logs/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_solver.log",
            python=RDIR + "/logs/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_python.log",
            memory=RDIR + "/logs/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_memory.log"
        threads: 4
        resources: mem_mb=config['solving']['mem']
        benchmark: RDIR + "/benchmarks/solve_network/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}"
        script: "scripts/solve_network.py"
