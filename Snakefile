from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()


configfile: "config.yaml"


SDIR = config["summary_dir"] + config["run"]
RDIR = config["results_dir"] + config["run"]
CDIR = config["costs_dir"]

CUTOUTS_PATH = (
    "cutouts/africa-2013-era5-tutorial.nc"
    if config["tutorial"]
    else "cutouts/africa-2013-era5.nc"
)


wildcard_constraints:
    ll="[a-z0-9\.]+",
    simpl="[a-zA-Z0-9]*|all",
    clusters="[0-9]+m?|all",
    opts="[-+a-zA-Z0-9]*",
    sopts="[-+a-zA-Z0-9\.\s]*",


subworkflow pypsaearth:
    workdir:
        "../pypsa-africa"
    snakefile:
        "../pypsa-africa/Snakefile"
    configfile:
        "./config.pypsa-earth.yaml"


rule prepare_sector_networks:
    input:
        expand(
            RDIR
            + "/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}.nc",
            **config["scenario"]
        ),


rule solve_all_networks:
    input:
        expand(
            RDIR
            + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}.nc",
            **config["scenario"]
        ),


rule prepare_sector_network:
    input:
        network=pypsaearth("networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"),
        costs=CDIR + "costs_{planning_horizons}.csv",
        h2_cavern="data/hydrogen_salt_cavern_potentials.csv",
        nodal_energy_totals="resources/nodal_energy_heat_totals_s{simpl}_{clusters}.csv",
        transport="resources/transport_s{simpl}_{clusters}.csv",
        avail_profile="resources/avail_profile_s{simpl}_{clusters}.csv",
        dsm_profile="resources/dsm_profile_s{simpl}_{clusters}.csv",
        nodal_transport_data="resources/nodal_transport_data_s{simpl}_{clusters}.csv",
        overrides="data/override_component_attrs",
        clustered_pop_layout="resources/pop_layout_elec_s{simpl}_{clusters}.csv",
        industrial_demand="resources/industrial_energy_demand_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
        airports="data/airports.csv",
        ports="data/ports.csv",
        heat_demand="resources/heat/heat_demand_s{simpl}_{clusters}.csv",
        ashp_cop="resources/heat/ashp_cop_s{simpl}_{clusters}.csv",
        gshp_cop="resources/heat/gshp_cop_s{simpl}_{clusters}.csv",
        solar_thermal="resources/heat/solar_thermal_s{simpl}_{clusters}.csv",
        district_heat_share="resources/heat/district_heat_share_s{simpl}_{clusters}.csv",
        industry_demands="data/industry_demand_locations.csv",
        biomass_potentials="data/temp_hard_coded/biomass_potentials_s_37.csv",
        biomass_transport_costs="data/temp_hard_coded/biomass_transport_costs.csv",
    output:
        RDIR
        + "/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}.nc",
    threads: 1
    resources:
        mem_mb=2000,
    benchmark:
        (
            RDIR
            + "/benchmarks/prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}"
        )
    script:
        "scripts/prepare_sector_network.py"


rule prepare_transport_data:
    input:
        network=pypsaearth("networks/elec_s{simpl}_{clusters}.nc"),
        energy_totals_name="resources/energy_totals.csv",
        traffic_data_KFZ="data/emobility/KFZ__count",
        traffic_data_Pkw="data/emobility/Pkw__count",
        transport_name="resources/transport_data.csv",
        clustered_pop_layout="resources/pop_layout_elec_s{simpl}_{clusters}.csv",
        temp_air_total="resources/temp_air_total_elec_s{simpl}_{clusters}.nc",
    output:
        nodal_energy_totals="resources/nodal_energy_totals_s{simpl}_{clusters}.csv",
        transport="resources/transport_s{simpl}_{clusters}.csv",
        avail_profile="resources/avail_profile_s{simpl}_{clusters}.csv",
        dsm_profile="resources/dsm_profile_s{simpl}_{clusters}.csv",
        nodal_transport_data="resources/nodal_transport_data_s{simpl}_{clusters}.csv",
    script:
        "scripts/prepare_transport_data.py"


rule build_cop_profiles:
    input:
        temp_soil_total="resources/temp_soil_total_elec_s{simpl}_{clusters}.nc",
        temp_soil_rural="resources/temp_soil_rural_elec_s{simpl}_{clusters}.nc",
        temp_soil_urban="resources/temp_soil_urban_elec_s{simpl}_{clusters}.nc",
        temp_air_total="resources/temp_air_total_elec_s{simpl}_{clusters}.nc",
        temp_air_rural="resources/temp_air_rural_elec_s{simpl}_{clusters}.nc",
        temp_air_urban="resources/temp_air_urban_elec_s{simpl}_{clusters}.nc",
    output:
        cop_soil_total="resources/cop_soil_total_elec_s{simpl}_{clusters}.nc",
        cop_soil_rural="resources/cop_soil_rural_elec_s{simpl}_{clusters}.nc",
        cop_soil_urban="resources/cop_soil_urban_elec_s{simpl}_{clusters}.nc",
        cop_air_total="resources/cop_air_total_elec_s{simpl}_{clusters}.nc",
        cop_air_rural="resources/cop_air_rural_elec_s{simpl}_{clusters}.nc",
        cop_air_urban="resources/cop_air_urban_elec_s{simpl}_{clusters}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        "benchmarks/build_cop_profiles/s{simpl}_{clusters}"
    script:
        "scripts/build_cop_profiles.py"


rule prepare_heat_data:
    input:
        network=pypsaearth("networks/elec_s{simpl}_{clusters}.nc"),
        energy_totals_name="resources/energy_totals.csv",
        clustered_pop_layout="resources/pop_layout_elec_s{simpl}_{clusters}.csv",
        temp_air_total="resources/temp_air_total_elec_s{simpl}_{clusters}.nc",
        cop_soil_total="resources/cop_soil_total_elec_s{simpl}_{clusters}.nc",
        cop_air_total="resources/cop_air_total_elec_s{simpl}_{clusters}.nc",
        solar_thermal_total="resources/solar_thermal_total_elec_s{simpl}_{clusters}.nc",
        heat_demand_total="resources/heat_demand_total_elec_s{simpl}_{clusters}.nc",
        heat_profile="data/heat_load_profile_BDEW.csv",
    output:
        nodal_energy_totals="resources/nodal_energy_heat_totals_s{simpl}_{clusters}.csv",
        heat_demand="resources/heat/heat_demand_s{simpl}_{clusters}.csv",
        ashp_cop="resources/heat/ashp_cop_s{simpl}_{clusters}.csv",
        gshp_cop="resources/heat/gshp_cop_s{simpl}_{clusters}.csv",
        solar_thermal="resources/heat/solar_thermal_s{simpl}_{clusters}.csv",
        district_heat_share="resources/heat/district_heat_share_s{simpl}_{clusters}.csv",
    script:
        "scripts/prepare_heat_data.py"


rule build_solar_thermal_profiles:
    input:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc",
        regions_onshore=pypsaearth(
            "resources/regions_onshore_elec_s{simpl}_{clusters}.geojson"
        ),
        cutout=pypsaearth(CUTOUTS_PATH),
    output:
        solar_thermal_total="resources/solar_thermal_total_elec_s{simpl}_{clusters}.nc",
        solar_thermal_urban="resources/solar_thermal_urban_elec_s{simpl}_{clusters}.nc",
        solar_thermal_rural="resources/solar_thermal_rural_elec_s{simpl}_{clusters}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        "benchmarks/build_solar_thermal_profiles/s{simpl}_{clusters}"
    script:
        "scripts/build_solar_thermal_profiles.py"


rule build_population_layouts:
    input:
        nuts3_shapes=pypsaearth("resources/gadm_shapes.geojson"),
        urban_percent="data/urban_percent.csv",
        cutout=pypsaearth(CUTOUTS_PATH),
    output:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc",
    resources:
        mem_mb=20000,
    benchmark:
        "benchmarks/build_population_layouts"
    threads: 8
    script:
        "scripts/build_population_layouts.py"


rule build_industrial_distribution_key:
    input:
        regions_onshore=pypsaearth(
            "resources/regions_onshore_elec_s{simpl}_{clusters}.geojson"
        ),
        clustered_pop_layout="resources/pop_layout_elec_s{simpl}_{clusters}.csv",
        industrial_database="data/morocco_cement_industry.csv",
    output:
        industrial_distribution_key="resources/industrial_distribution_key_elec_s{simpl}_{clusters}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/build_industrial_distribution_key/s{simpl}_{clusters}"
    script:
        "scripts/build_industrial_distribution_key.py"


rule build_industrial_energy_demand_per_node:
    input:
        industry_sector_ratios="data/industry_sector_ratios.csv",
        industrial_energy_demand_per_node_today="resources/industrial_energy_demand_today_elec_s{simpl}_{clusters}.csv",
        industrial_distribution_key="resources/industrial_distribution_key_elec_s{simpl}_{clusters}.csv",
        industrial_production_per_country_tomorrow="data/industrial_production_per_country_tomorrow_{planning_horizons}.csv",
    output:
        industrial_energy_demand_per_node="resources/industrial_energy_demand_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/build_industrial_energy_demand_per_node/s{simpl}_{clusters}_{planning_horizons}"
    script:
        "scripts/build_industrial_energy_demand_per_node.py"


rule build_industrial_energy_demand_per_node_today:
    input:
        industrial_distribution_key="resources/industrial_distribution_key_elec_s{simpl}_{clusters}.csv",
        industrial_energy_demand_per_country_today="data/industrial_energy_demand_per_country_today.csv",
    output:
        industrial_energy_demand_per_node_today="resources/industrial_energy_demand_today_elec_s{simpl}_{clusters}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/build_industrial_energy_demand_per_node_today/s{simpl}_{clusters}"
    script:
        "scripts/build_industrial_energy_demand_per_node_today.py"


rule move_hardcoded_files_temp:
    input:
        "data/temp_hard_coded/energy_totals.csv",
        "data/temp_hard_coded/transport_data.csv",
    output:
        "resources/energy_totals.csv",
        "resources/transport_data.csv",
    shell:
        "cp -a data/temp_hard_coded/. resources"


rule build_clustered_population_layouts:
    input:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc",
        regions_onshore=pypsaearth(
            "resources/regions_onshore_elec_s{simpl}_{clusters}.geojson"
        ),
        cutout=pypsaearth(CUTOUTS_PATH),
    output:
        clustered_pop_layout="resources/pop_layout_elec_s{simpl}_{clusters}.csv",
    resources:
        mem_mb=10000,
    benchmark:
        "benchmarks/build_clustered_population_layouts/s{simpl}_{clusters}"
    script:
        "scripts/build_clustered_population_layouts.py"


rule build_heat_demand:
    input:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc",
        regions_onshore=pypsaearth(
            "resources/regions_onshore_elec_s{simpl}_{clusters}.geojson"
        ),
        cutout=pypsaearth(CUTOUTS_PATH),
    output:
        heat_demand_urban="resources/heat_demand_urban_elec_s{simpl}_{clusters}.nc",
        heat_demand_rural="resources/heat_demand_rural_elec_s{simpl}_{clusters}.nc",
        heat_demand_total="resources/heat_demand_total_elec_s{simpl}_{clusters}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        "benchmarks/build_heat_demand/s{simpl}_{clusters}"
    script:
        "scripts/build_heat_demand.py"


rule build_temperature_profiles:
    input:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc",
        regions_onshore=pypsaearth(
            "resources/regions_onshore_elec_s{simpl}_{clusters}.geojson"
        ),
        cutout=pypsaearth(CUTOUTS_PATH),
    output:
        temp_soil_total="resources/temp_soil_total_elec_s{simpl}_{clusters}.nc",
        temp_soil_rural="resources/temp_soil_rural_elec_s{simpl}_{clusters}.nc",
        temp_soil_urban="resources/temp_soil_urban_elec_s{simpl}_{clusters}.nc",
        temp_air_total="resources/temp_air_total_elec_s{simpl}_{clusters}.nc",
        temp_air_rural="resources/temp_air_rural_elec_s{simpl}_{clusters}.nc",
        temp_air_urban="resources/temp_air_urban_elec_s{simpl}_{clusters}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        "benchmarks/build_temperature_profiles/s{simpl}_{clusters}"
    script:
        "scripts/build_temperature_profiles.py"


rule copy_config:
    output:
        SDIR + "/configs/config.yaml",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        SDIR + "/benchmarks/copy_config"
    script:
        "scripts/copy_config.py"


rule solve_network:
    input:
        overrides="data/override_component_attrs",
        network=RDIR
        + "/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}.nc",
        costs=CDIR + "costs_{planning_horizons}.csv",
    output:
        RDIR
        + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}.nc",
    shadow:
        "shallow"
    log:
        solver=RDIR
        + "/logs/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_solver.log",
        python=RDIR
        + "/logs/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_python.log",
        memory=RDIR
        + "/logs/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_memory.log",
    threads: 4
    resources:
        mem_mb=config["solving"]["mem"],
    benchmark:
        (
            RDIR
            + "/benchmarks/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}"
        )
    script:
        "scripts/solve_network.py"


rule make_summary:
    input:
        overrides="data/override_component_attrs",
        networks=expand(
            RDIR
            + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}.nc",
            **config["scenario"]
        ),
        costs=CDIR + "costs_{}.csv".format(config["scenario"]["planning_horizons"][0]),
        plots=expand(
            RDIR
            + "/maps/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}-costs-all_{planning_horizons}.pdf",
            **config["scenario"]
        ),
    output:
        nodal_costs=SDIR + "/csvs/nodal_costs.csv",
        nodal_capacities=SDIR + "/csvs/nodal_capacities.csv",
        nodal_cfs=SDIR + "/csvs/nodal_cfs.csv",
        cfs=SDIR + "/csvs/cfs.csv",
        costs=SDIR + "/csvs/costs.csv",
        capacities=SDIR + "/csvs/capacities.csv",
        curtailment=SDIR + "/csvs/curtailment.csv",
        energy=SDIR + "/csvs/energy.csv",
        supply=SDIR + "/csvs/supply.csv",
        supply_energy=SDIR + "/csvs/supply_energy.csv",
        prices=SDIR + "/csvs/prices.csv",
        weighted_prices=SDIR + "/csvs/weighted_prices.csv",
        market_values=SDIR + "/csvs/market_values.csv",
        price_statistics=SDIR + "/csvs/price_statistics.csv",
        metrics=SDIR + "/csvs/metrics.csv",
    threads: 2
    resources:
        mem_mb=10000,
    benchmark:
        SDIR + "/benchmarks/make_summary"
    script:
        "scripts/make_summary.py"


rule plot_network:
    input:
        overrides="data/override_component_attrs",
        network=RDIR
        + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}.nc",
    output:
        map=RDIR
        + "/maps/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}-costs-all_{planning_horizons}.pdf",
    threads: 2
    resources:
        mem_mb=10000,
    benchmark:
        (
            RDIR
            + "/benchmarks/plot_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}"
        )
    script:
        "scripts/plot_network.py"


rule plot_summary:
    input:
        costs=SDIR + "/csvs/costs.csv",
        energy=SDIR + "/csvs/energy.csv",
        balances=SDIR + "/csvs/supply_energy.csv",
    output:
        costs=SDIR + "/graphs/costs.pdf",
        energy=SDIR + "/graphs/energy.pdf",
        balances=SDIR + "/graphs/balances-energy.pdf",
    threads: 2
    resources:
        mem_mb=10000,
    benchmark:
        SDIR + "/benchmarks/plot_summary"
    script:
        "scripts/plot_summary.py"


rule prepare_db:
    input:
        network=RDIR
        + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}.nc",
    output:
        db=RDIR
        + "/summaries/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}-costs-all_{planning_horizons}.csv",
    threads: 2
    resources:
        mem_mb=10000,
    benchmark:
        (
            RDIR
            + "/benchmarks/prepare_db/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}"
        )
    script:
        "scripts/prepare_db.py"
