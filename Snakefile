
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

configfile: "config.yaml"

SDIR = config['summary_dir'] + config['run']
RDIR = config['results_dir'] + config['run']
CDIR = config['costs_dir']

wildcard_constraints:
    lv="[a-z0-9\.]+",
    simpl="[a-zA-Z0-9]*|all",
    clusters="[0-9]+m?|all",
    opts="[-+a-zA-Z0-9]*",
    sector_opts="[-+a-zA-Z0-9\.\s]*"

rule prepare_sector_networks:
    input:
        expand(RDIR + "/prenetworks/elec_s{simpl}_{clusters}_{planning_horizons}.nc",
               **config['scenario'])


rule solve_all_networks:
    input:
        expand(RDIR + "/postnetworks/elec_s{simpl}_{clusters}_{planning_horizons}.nc",
               **config['scenario'])


rule prepare_sector_network:
    input:
        network='networks/elec_s{simpl}_{clusters}.nc',
        costs=CDIR + "costs_{planning_horizons}.csv",
        h2_cavern="data/hydrogen_salt_cavern_potentials.csv",
        nodal_energy_totals='resources/nodal_energy_totals.csv',
        transport='resources/transport.csv',
        avail_profile='resources/avail_profile.csv',
        dsm_profile='resources/dsm_profile.csv',
        nodal_transport_data='resources/nodal_transport_data.csv',
        overrides="data/override_component_attrs",

    output: RDIR + '/prenetworks/elec_s{simpl}_{clusters}_{planning_horizons}.nc'
    threads: 1
    resources: mem_mb=2000
    benchmark: RDIR + "/benchmarks/prepare_network/elec_s{simpl}_{clusters}_{planning_horizons}"
    script: "scripts/prepare_sector_network.py"

rule prepare_transport_data:
    input:
        network='networks/elec_s_4.nc',  # hardcoded wildcards temporary
        # network='networks/elec_s{simpl}_{clusters}.nc',
        energy_totals_name='resources/energy_totals.csv',
        traffic_data_KFZ = "data/emobility/KFZ__count",
        traffic_data_Pkw = "data/emobility/Pkw__count",
        transport_name='resources/transport_data.csv',
        clustered_pop_layout="resources/pop_layout_elec_s_4.csv",  # hardcoded wildcards temporary
        # clustered_pop_layout="resources/pop_layout_elec_s{simpl}_{clusters}.csv",
        # This is probably still dummy data, investigate and use real data TODO
        temp_air_total="resources/temp_air_total_elec_s_4.nc",  # hardcoded wildcards temporary
        # temp_air_total="resources/temp_air_total_elec_s{simpl}_{clusters}.nc",

    output: 
        nodal_energy_totals='resources/nodal_energy_totals.csv',
        transport='resources/transport.csv',
        avail_profile='resources/avail_profile.csv',
        dsm_profile='resources/dsm_profile.csv',
        nodal_transport_data='resources/nodal_transport_data.csv',
        # dummy_wildcard="resources/dummy{simpl}_{clusters}.nc"

    script: "scripts/prepare_transport_data.py"

rule calculate_dummy_pop_layout:
    input:
        network='networks/elec_s{simpl}_37.nc',  # fixed number of clusters
        # Get pop layouts from Europe (update to Morocco/Africa layout)
        clustered_pop_layout="resources/pop_layout_elec_s{simpl}_37.csv",  # fixed number of clusters
        #simplified_pop_layout="resources/pop_layout_elec_s{simpl}.csv",

    output: clustered_pop_layout_dummy="resources/pop_layout_elec_s{simpl}_37.csv",  # fixed number of clusters
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

rule build_clustered_population_layouts:
    input:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc",
        regions_onshore="resources/regions_onshore_elec_s{simpl}_{clusters}.geojson"
    output:
        clustered_pop_layout="resources/pop_layout_elec_s{simpl}_{clusters}.csv"
    resources: mem_mb=10000
    benchmark: "benchmarks/build_clustered_population_layouts/s{simpl}_{clusters}"
    script: "scripts/build_clustered_population_layouts.py"
    
    
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
        network=RDIR + "/prenetworks/elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        costs=CDIR + "costs_{planning_horizons}.csv",
        
    output: RDIR + "/postnetworks/elec_s{simpl}_{clusters}_{planning_horizons}.nc"
    shadow: "shallow"
    log:
        solver=RDIR + "/logs/elec_s{simpl}_{clusters}_{planning_horizons}_solver.log",
        python=RDIR + "/logs/elec_s{simpl}_{clusters}_{planning_horizons}_python.log",
        memory=RDIR + "/logs/elec_s{simpl}_{clusters}_{planning_horizons}_memory.log",
    threads: 4
    resources: mem_mb=config['solving']['mem']
    benchmark: RDIR + "/benchmarks/solve_network/elec_s{simpl}_{clusters}_{planning_horizons}"
    script: "scripts/solve_network.py"

rule make_summary:
    input:
        overrides="data/override_component_attrs",
        networks=expand(
            RDIR + "/postnetworks/elec_s{simpl}_{clusters}_{planning_horizons}.nc",
            **config['scenario']
        ),
        costs=CDIR + "costs_{}.csv".format(config['scenario']['planning_horizons'][0]),
        plots=expand(
            RDIR + "/maps/elec_s{simpl}_{clusters}-costs-all_{planning_horizons}.pdf",
            **config['scenario']
        )
    output:
        nodal_costs=SDIR + '/csvs/nodal_costs.csv',
        nodal_capacities=SDIR + '/csvs/nodal_capacities.csv',
        nodal_cfs=SDIR + '/csvs/nodal_cfs.csv',
        cfs=SDIR + '/csvs/cfs.csv',
        costs=SDIR + '/csvs/costs.csv',
        capacities=SDIR + '/csvs/capacities.csv',
        curtailment=SDIR + '/csvs/curtailment.csv',
        energy=SDIR + '/csvs/energy.csv',
        supply=SDIR + '/csvs/supply.csv',
        supply_energy=SDIR + '/csvs/supply_energy.csv',
        prices=SDIR + '/csvs/prices.csv',
        weighted_prices=SDIR + '/csvs/weighted_prices.csv',
        market_values=SDIR + '/csvs/market_values.csv',
        price_statistics=SDIR + '/csvs/price_statistics.csv',
        metrics=SDIR + '/csvs/metrics.csv'
    threads: 2
    resources: mem_mb=10000
    benchmark: SDIR + "/benchmarks/make_summary"
    script: "scripts/make_summary.py"

rule plot_network:
    input:
        overrides="data/override_component_attrs",
        network=RDIR + "/postnetworks/elec_s{simpl}_{clusters}_{planning_horizons}.nc"
    output:
        map=RDIR + "/maps/elec_s{simpl}_{clusters}_-costs-all_{planning_horizons}.pdf",
        today=RDIR + "/maps/elec_s{simpl}_{clusters}_{planning_horizons}-today.pdf"
    threads: 2
    resources: mem_mb=10000
    benchmark: RDIR + "/benchmarks/plot_network/elec_s{simpl}_{clusters}_{planning_horizons}"
    script: "scripts/plot_network_plot_network_pypsa_earth.py"

rule plot_summary:
    input:
        costs=SDIR + '/csvs/costs.csv',
        energy=SDIR + '/csvs/energy.csv',
        balances=SDIR + '/csvs/supply_energy.csv'
    output:
        costs=SDIR + '/graphs/costs.pdf',
        energy=SDIR + '/graphs/energy.pdf',
        balances=SDIR + '/graphs/balances-energy.pdf'
    threads: 2
    resources: mem_mb=10000
    benchmark: SDIR + "/benchmarks/plot_summary"
    script: "scripts/plot_summary.py"
    