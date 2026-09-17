# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later


solar_rooftop_config = config["sector"]["solar_rooftop"]
if isinstance(solar_rooftop_config, dict):
    solar_rooftop_enable = (
        solar_rooftop_config["enable"] and solar_rooftop_config["use_building_size"]
    )
    solar_rooftop_params = {
        "solar_rooftop_enable": solar_rooftop_enable,
        "install_ratio": solar_rooftop_config["install_ratio"],
        "tolerance": solar_rooftop_config["tolerance"],
    }
else:
    solar_rooftop_params = {}
    solar_rooftop_enable = config["sector"]["solar_rooftop"]


rule cluster_global_buildings:
    params:
        **solar_rooftop_params,
        crs=config["crs"],
    input:
        country_buildings="data/global_buildings/{country}_global_buildings_raw.parquet",
        regions_onshore=rules.cluster_network.output.regions_onshore,
    output:
        solar_rooftop_layout=branch(
            solar_rooftop_enable,
            "resources/"
            + RDIR
            + "solar_rooftop/solar_rooftop_layout_elec_s{simpl}_{clusters}_{country}.csv",
        ),
    script:
        scripts("cluster_global_buildings.py")


rule prepare_ports:
    params:
        custom_export=config["custom_data"]["export_ports"],
    output:
        ports="resources/" + SECDIR + "ports.csv",
        export_ports="resources/" + SECDIR + "export_ports.csv",
    script:
        scripts("prepare_ports.py")


rule prepare_airports:
    params:
        airport_sizing_factor=config["sector"]["airport_sizing_factor"],
        airport_custom_data=config["custom_data"]["airports"],
    output:
        airports="resources/" + SECDIR + "airports.csv",
    script:
        scripts("prepare_airports.py")


rule prepare_urban_percent:
    output:
        urban_percent="resources/" + SECDIR + "urban_percent.csv",
    script:
        scripts("prepare_urban_percent.py")


rule prepare_transport_data_input:
    output:
        transport_data_input="resources/" + SECDIR + "transport_data.csv",
    script:
        scripts("prepare_transport_data_input.py")


if (
    not config["custom_data"]["h2_underground"]
    and config["sector"]["hydrogen"]["underground_storage"]["enabled"]
):

    rule build_salt_cavern_potentials:
        input:
            copernicus="data/copernicus/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
            regions_onshore=rules.cluster_network.output.regions_onshore,
            regions_offshore=rules.cluster_network.output.regions_offshore,
            potash_shp=rules.retrieve_potash_data.output.potash_files,
        output:
            h2_cavern="resources/"
            + RDIR
            + "salt_cavern_potentials_s{simpl}_{clusters}.csv",
        params:
            crs=config["crs"],
            underground_storage=config["sector"]["hydrogen"]["underground_storage"],
        threads: 1
        resources:
            mem_mb=2000,
        script:
            scripts("build_salt_cavern_potentials.py")


if not config["custom_data"]["gas_network"]:

    rule prepare_gas_network:
        params:
            gas_config=config["sector"]["gas"],
            alternative_clustering=config["clustering"]["alternative_clustering"],
            custom_gas_network=config["custom_data"]["gas_network"],
        input:
            regions_onshore=rules.cluster_network.output.regions_onshore,
        output:
            clustered_gas_network="resources/"
            + SECDIR
            + "gas_networks/gas_network_elec_s{simpl}_{clusters}.csv",
            # TODO: Should be a own snakemake rule
            # gas_network_fig_1="resources/gas_networks/existing_gas_pipelines_{simpl}_{clusters}.png",
            # gas_network_fig_2="resources/gas_networks/clustered_gas_pipelines_{simpl}_{clusters}.png",
        script:
            scripts("prepare_gas_network.py")


rule build_population_layouts:
    params:
        planning_horizons=config["scenario"]["planning_horizons"][0],
    input:
        nuts3_shapes=rules.build_shapes.output.gadm_shapes,
        urban_percent=rules.prepare_urban_percent.output.urban_percent,
        cutout="cutouts/"
        + CDIR
        + [c["cutout"] for _, c in config["renewable"].items()][0]
        + ".nc",
        # default to first cutout found
    output:
        pop_layout_total="resources/"
        + SECDIR
        + "population_shares/pop_layout_total_{planning_horizons}.nc",
        pop_layout_urban="resources/"
        + SECDIR
        + "population_shares/pop_layout_urban_{planning_horizons}.nc",
        pop_layout_rural="resources/"
        + SECDIR
        + "population_shares/pop_layout_rural_{planning_horizons}.nc",
        gdp_layout="resources/"
        + SECDIR
        + "gdp_shares/gdp_layout_{planning_horizons}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        ("benchmarks/" + SECDIR + "build_population_layouts_{planning_horizons}")
    threads: 8
    script:
        scripts("build_population_layouts.py")


rule build_clustered_population_layouts:
    input:
        pop_layout_total=rules.build_population_layouts.output.pop_layout_total,
        pop_layout_urban=rules.build_population_layouts.output.pop_layout_urban,
        pop_layout_rural=rules.build_population_layouts.output.pop_layout_rural,
        gdp_layout=rules.build_population_layouts.output.gdp_layout,
        regions_onshore=rules.cluster_network.output.regions_onshore,
        cutout="cutouts/"
        + CDIR
        + [c["cutout"] for _, c in config["renewable"].items()][0]
        + ".nc",
        # default to first cutout found
    output:
        clustered_pop_layout="resources/"
        + SECDIR
        + "population_shares/pop_layout_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
        clustered_gdp_layout="resources/"
        + SECDIR
        + "gdp_shares/gdp_layout_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
    resources:
        mem_mb=10000,
    benchmark:
        (
            "benchmarks/"
            + SECDIR
            + "build_clustered_population_layouts/s{simpl}_{clusters}_{planning_horizons}"
        )
    script:
        scripts("build_clustered_population_layouts.py")


rule build_base_energy_totals:
    params:
        space_heat_share=config["sector"]["space_heat_share"],
        update_data=config["demand_data"]["update_data"],
        base_year=config["demand_data"]["base_year"],
        countries=config["countries"],
        shift_coal_to_elec=config["sector"]["coal"]["shift_to_elec"],
    input:
        unsd_paths="data/demand/unsd/paths/Energy_Statistics_Database.xlsx",
    output:
        energy_totals_base="resources/" + SECDIR + "energy_totals_base.csv",
        unsd_export_path=directory("data/demand/unsd/data/"),
    script:
        scripts("build_base_energy_totals.py")


rule prepare_energy_totals:
    params:
        countries=config["countries"],
        base_year=config["demand_data"]["base_year"],
        sector_options=config["sector"],
        demand_scenario=config["demand_data"]["scenario"],
    input:
        unsd_paths=rules.build_base_energy_totals.output.energy_totals_base,
        efficiency_gains_cagr="data/demand/efficiency_gains_cagr.csv",
        growth_factors_cagr="data/demand/growth_factors_cagr.csv",
        district_heating="data/demand/district_heating.csv",
        fuel_shares="data/demand/fuel_shares.csv",
    output:
        energy_totals="resources/" + SECDIR + "energy_totals_{planning_horizons}.csv",
    script:
        scripts("prepare_energy_totals.py")


rule build_heat_demand:
    params:
        snapshots=config["snapshots"],
    input:
        pop_layout_total=rules.build_population_layouts.output.pop_layout_total,
        pop_layout_urban=rules.build_population_layouts.output.pop_layout_urban,
        pop_layout_rural=rules.build_population_layouts.output.pop_layout_rural,
        regions_onshore=rules.cluster_network.output.regions_onshore,
        cutout="cutouts/"
        + CDIR
        + [c["cutout"] for _, c in config["renewable"].items()][0]
        + ".nc",
        # default to first cutout found
    output:
        heat_demand_urban="resources/"
        + SECDIR
        + "demand/heat/heat_demand_urban_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        heat_demand_rural="resources/"
        + SECDIR
        + "demand/heat/heat_demand_rural_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        heat_demand_total="resources/"
        + SECDIR
        + "demand/heat/heat_demand_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        (
            "benchmarks/"
            + SECDIR
            + "build_heat_demand/s{simpl}_{clusters}_{planning_horizons}"
        )
    script:
        scripts("build_heat_demand.py")


rule build_temperature_profiles:
    params:
        snapshots=config["snapshots"],
    input:
        pop_layout_total=rules.build_population_layouts.output.pop_layout_total,
        pop_layout_urban=rules.build_population_layouts.output.pop_layout_urban,
        pop_layout_rural=rules.build_population_layouts.output.pop_layout_rural,
        regions_onshore=rules.cluster_network.output.regions_onshore,
        cutout="cutouts/"
        + CDIR
        + [c["cutout"] for _, c in config["renewable"].items()][0]
        + ".nc",
        # default to first cutout found
    output:
        temp_soil_total="resources/"
        + SECDIR
        + "temperatures/temp_soil_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        temp_soil_rural="resources/"
        + SECDIR
        + "temperatures/temp_soil_rural_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        temp_soil_urban="resources/"
        + SECDIR
        + "temperatures/temp_soil_urban_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        temp_air_total="resources/"
        + SECDIR
        + "temperatures/temp_air_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        temp_air_rural="resources/"
        + SECDIR
        + "temperatures/temp_air_rural_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        temp_air_urban="resources/"
        + SECDIR
        + "temperatures/temp_air_urban_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        (
            "benchmarks/"
            + SECDIR
            + "build_temperature_profiles/s{simpl}_{clusters}_{planning_horizons}"
        )
    script:
        scripts("build_temperature_profiles.py")


rule build_cop_profiles:
    params:
        heat_pump_sink_T=config["sector"]["heat_pump_sink_T"],
    input:
        temp_soil_total=rules.build_temperature_profiles.output.temp_soil_total,
        temp_soil_rural=rules.build_temperature_profiles.output.temp_soil_rural,
        temp_soil_urban=rules.build_temperature_profiles.output.temp_soil_urban,
        temp_air_total=rules.build_temperature_profiles.output.temp_air_total,
        temp_air_rural=rules.build_temperature_profiles.output.temp_air_rural,
        temp_air_urban=rules.build_temperature_profiles.output.temp_air_urban,
    output:
        cop_soil_total="resources/"
        + SECDIR
        + "cops/cop_soil_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        cop_soil_rural="resources/"
        + SECDIR
        + "cops/cop_soil_rural_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        cop_soil_urban="resources/"
        + SECDIR
        + "cops/cop_soil_urban_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        cop_air_total="resources/"
        + SECDIR
        + "cops/cop_air_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        cop_air_rural="resources/"
        + SECDIR
        + "cops/cop_air_rural_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        cop_air_urban="resources/"
        + SECDIR
        + "cops/cop_air_urban_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        (
            "benchmarks/"
            + SECDIR
            + "build_cop_profiles/s{simpl}_{clusters}_{planning_horizons}"
        )
    script:
        scripts("build_cop_profiles.py")


rule build_solar_thermal_profiles:
    params:
        solar_thermal_config=config["sector"]["solar_thermal_collector"],
        snapshots=config["snapshots"],
    input:
        pop_layout_total=rules.build_population_layouts.output.pop_layout_total,
        pop_layout_urban=rules.build_population_layouts.output.pop_layout_urban,
        pop_layout_rural=rules.build_population_layouts.output.pop_layout_rural,
        regions_onshore=rules.cluster_network.output.regions_onshore,
        cutout="cutouts/"
        + CDIR
        + [c["cutout"] for _, c in config["renewable"].items()][0]
        + ".nc",
        # default to first cutout found
    output:
        solar_thermal_total="resources/"
        + SECDIR
        + "demand/heat/solar_thermal_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        solar_thermal_urban="resources/"
        + SECDIR
        + "demand/heat/solar_thermal_urban_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        solar_thermal_rural="resources/"
        + SECDIR
        + "demand/heat/solar_thermal_rural_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        (
            "benchmarks/"
            + SECDIR
            + "build_solar_thermal_profiles/s{simpl}_{clusters}_{planning_horizons}"
        )
    script:
        scripts("build_solar_thermal_profiles.py")


rule prepare_heat_data:
    input:
        network=clustered_network,
        energy_totals_name=rules.prepare_energy_totals.output.energy_totals,
        clustered_pop_layout=rules.build_clustered_population_layouts.output.clustered_pop_layout,
        temp_air_total=rules.build_temperature_profiles.output.temp_air_total,
        cop_soil_total=rules.build_cop_profiles.output.cop_soil_total,
        cop_air_total=rules.build_cop_profiles.output.cop_air_total,
        solar_thermal_total=rules.build_solar_thermal_profiles.output.solar_thermal_total,
        heat_demand_total=rules.build_heat_demand.output.heat_demand_total,
        heat_profile="data/heat_load_profile_BDEW.csv",
    output:
        nodal_energy_totals="resources/"
        + SECDIR
        + "demand/heat/nodal_energy_heat_totals_s{simpl}_{clusters}_{planning_horizons}.csv",
        heat_demand="resources/"
        + SECDIR
        + "demand/heat/heat_demand_s{simpl}_{clusters}_{planning_horizons}.csv",
        ashp_cop="resources/"
        + SECDIR
        + "demand/heat/ashp_cop_s{simpl}_{clusters}_{planning_horizons}.csv",
        gshp_cop="resources/"
        + SECDIR
        + "demand/heat/gshp_cop_s{simpl}_{clusters}_{planning_horizons}.csv",
        solar_thermal="resources/"
        + SECDIR
        + "demand/heat/solar_thermal_s{simpl}_{clusters}_{planning_horizons}.csv",
        district_heat_share="resources/"
        + SECDIR
        + "demand/heat/district_heat_share_s{simpl}_{clusters}_{planning_horizons}.csv",
    script:
        scripts("prepare_heat_data.py")


rule prepare_transport_data:
    input:
        network=clustered_network,
        energy_totals_name=rules.prepare_energy_totals.output.energy_totals,
        traffic_data_KFZ="data/emobility/KFZ__count",
        traffic_data_Pkw="data/emobility/Pkw__count",
        transport_name=rules.prepare_transport_data_input.output.transport_data_input,
        clustered_pop_layout=rules.build_clustered_population_layouts.output.clustered_pop_layout,
        temp_air_total=rules.build_temperature_profiles.output.temp_air_total,
    output:
        # nodal_energy_totals="resources/nodal_energy_totals_s{simpl}_{clusters}.csv",
        transport="resources/"
        + SECDIR
        + "demand/transport_s{simpl}_{clusters}_{planning_horizons}.csv",
        avail_profile="resources/"
        + SECDIR
        + "pattern_profiles/avail_profile_s{simpl}_{clusters}_{planning_horizons}.csv",
        dsm_profile="resources/"
        + SECDIR
        + "pattern_profiles/dsm_profile_s{simpl}_{clusters}_{planning_horizons}.csv",
        nodal_transport_data="resources/"
        + SECDIR
        + "demand/nodal_transport_data_s{simpl}_{clusters}_{planning_horizons}.csv",
    script:
        scripts("prepare_transport_data.py")


rule move_hardcoded_files_temp:
    input:
        "data/temp_hard_coded/energy_totals.csv",
    output:
        "resources/" + SECDIR + "energy_totals.csv",
    shell:
        "cp -a data/temp_hard_coded/. resources"


rule build_ammonia_production:
    input:
        ammonia_plants="data/industry/ammonia_plants.csv",
        us_cities=rules.retrieve_us_cities_dataset.output.us_cities,
        usgs_ammonia_dataset=rules.retrieve_ammonia_dataset.output.usgs_ammonia_dataset,
    output:
        ammonia_production="resources/ammonia_production.csv",
        ammonia_plants="resources/ammonia_plants.csv",
    threads: 1
    resources:
        mem_mb=1000,
    log:
        RESDIR + "logs/build_ammonia_production.log",
    benchmark:
        RESDIR + "benchmarks/build_ammonia_production"
    script:
        scripts("build_ammonia_production.py")


rule build_industrial_database:
    input:
        ammonia_plants=rules.build_ammonia_production.output.ammonia_plants,
    output:
        industrial_database="resources/industrial_database.csv",
    script:
        scripts("build_industrial_database.py")


rule build_industrial_distribution_key:  #default data
    params:
        countries=config["countries"],
        gadm_layer_id=config["build_shape_options"]["gadm_layer_id"],
        alternative_clustering=config["clustering"]["alternative_clustering"],
        industry_database=config["custom_data"]["industry_database"],
    input:
        regions_onshore=rules.cluster_network.output.regions_onshore,
        clustered_pop_layout=rules.build_clustered_population_layouts.output.clustered_pop_layout,
        clustered_gdp_layout=rules.build_clustered_population_layouts.output.clustered_gdp_layout,
        industrial_database=rules.build_industrial_database.output.industrial_database,
        shapes_path=rules.cluster_network.output.regions_onshore,
    output:
        industrial_distribution_key="resources/"
        + SECDIR
        + "demand/industrial_distribution_key_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        (
            "benchmarks/"
            + RDIR
            + "build_industrial_distribution_key_elec_s{simpl}_{clusters}_{planning_horizons}"
        )
    script:
        scripts("build_industrial_distribution_key.py")


rule build_base_industry_totals:  #default data
    params:
        base_year=config["demand_data"]["base_year"],
        countries=config["countries"],
        other_industries=config["demand_data"]["other_industries"],
        demand_scenario=config["demand_data"]["scenario"],
    input:
        #os.path.dirname(snakemake.input["transactions_path"]) + "/demand/unsd/data/"
        #industrial_production_per_country="data/industrial_production_per_country.csv",
        unsd_export_path="data/demand/unsd/data/",
        energy_totals_base=rules.build_base_energy_totals.output.energy_totals_base,
        transactions_path="data/unsd_transactions.csv",
    output:
        base_industry_totals="resources/"
        + SECDIR
        + "demand/base_industry_totals_{planning_horizons}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        ("benchmarks/" + SECDIR + "build_base_industry_totals_{planning_horizons}")
    script:
        scripts("build_base_industry_totals.py")


rule build_industry_demand:  #default data
    params:
        countries=config["countries"],
        industry_demand=config["custom_data"]["industry_demand"],
        demand_scenario=config["demand_data"]["scenario"],
        base_year=config["demand_data"]["base_year"],
        industry_util_factor=config["sector"]["industry_util_factor"],
        aluminium_year=config["demand_data"]["aluminium_year"],
        ammonia_enable=config["sector"]["ammonia"]["enable"],
        ammonia_gas_mwh_per_t=config["sector"]["ammonia"]["gas_MWh_per_tNH3"],
        ammonia_elec_mwh_per_t=config["sector"]["ammonia"]["elec_MWh_per_tNH3"],
        ammonia_year=config["sector"]["ammonia"]["production_year"],
    input:
        industrial_distribution_key=rules.build_industrial_distribution_key.output.industrial_distribution_key,
        #industrial_production_per_country_tomorrow="resources/demand/industrial_production_per_country_tomorrow_{planning_horizons}.csv",
        #industrial_production_per_country="data/industrial_production_per_country.csv",
        base_industry_totals=rules.build_base_industry_totals.output.base_industry_totals,
        industrial_database=rules.build_industrial_database.output.industrial_database,
        ammonia_production=rules.build_ammonia_production.output.ammonia_production,
        costs=rules.process_cost_data.output.costs.format(
            year="{planning_horizons}", scope="sec"
        ),
        industry_growth_cagr="data/demand/industry_growth_cagr.csv",
    output:
        industrial_energy_demand_per_node="resources/"
        + SECDIR
        + "demand/industrial_energy_demand_per_node_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        (
            "benchmarks/"
            + SECDIR
            + "industrial_energy_demand_per_node_elec_s{simpl}_{clusters}_{planning_horizons}.csv"
        )
    script:
        scripts("build_industry_demand.py")


rule build_existing_heating_distribution:
    params:
        baseyear=config["scenario"]["planning_horizons"][0],
        sector=config["sector"],
        existing_capacities=config["existing_capacities"],
    input:
        existing_heating="data/existing_infrastructure/existing_heating_raw.csv",
        clustered_pop_layout=rules.build_clustered_population_layouts.output.clustered_pop_layout,
        clustered_pop_energy_layout=rules.prepare_heat_data.output.nodal_energy_totals,
        #"resources/population_shares/pop_weighted_energy_totals_s{simpl}_{clusters}.csv",
        district_heat_share=rules.prepare_heat_data.output.district_heat_share,
    output:
        existing_heating_distribution="resources/"
        + SECDIR
        + "heating/existing_heating_distribution_s{simpl}_{clusters}_{planning_horizons}.csv",
    threads: 1
    resources:
        mem_mb=2000,
    log:
        RESDIR
        + "logs/build_existing_heating_distribution_s{simpl}_{clusters}_{planning_horizons}.log",
    benchmark:
        RESDIR
        +"benchmarks/build_existing_heating_distribution/s{simpl}_{clusters}_{planning_horizons}"
    script:
        scripts("build_existing_heating_distribution.py")


sector_enable = config["sector"]["enable"]

TRANSPORT = {
    "transport": rules.prepare_transport_data.output.transport,
    "avail_profile": rules.prepare_transport_data.output.avail_profile,
    "dsm_profile": rules.prepare_transport_data.output.dsm_profile,
    "nodal_transport_data": rules.prepare_transport_data.output.nodal_transport_data,
}

HEAT = {
    "heat_demand": rules.prepare_heat_data.output.heat_demand,
    "ashp_cop": rules.prepare_heat_data.output.ashp_cop,
    "gshp_cop": rules.prepare_heat_data.output.gshp_cop,
    "solar_thermal": rules.prepare_heat_data.output.solar_thermal,
    "district_heat_share": rules.prepare_heat_data.output.district_heat_share,
}


rule prepare_sector_network:
    params:
        electricity=config["electricity"],
        h2_underground=config["custom_data"]["h2_underground"],
        countries=config["countries"],
        gadm_layer_id=config["build_shape_options"]["gadm_layer_id"],
        alternative_clustering=config["clustering"]["alternative_clustering"],
        h2_policy=config["policy_config"]["hydrogen"],
        sector_options=config["sector"],
        foresight=config["foresight"],
        water_costs=config["custom_data"]["water_costs"],
        co2=config["co2"],
        demand_scenario=config["demand_data"]["scenario"],
    input:
        **branch(sector_enable["land_transport"], TRANSPORT),
        **branch(sector_enable["heat"], HEAT),
        **branch(
            config["custom_data"]["h2_underground"]
            or config["sector"]["hydrogen"]["underground_storage"]["enabled"],
            {
                "h2_cavern": branch(
                    config["custom_data"]["h2_underground"],
                    "data/hydrogen_salt_cavern_potentials.csv",
                    f"resources/{RDIR}salt_cavern_potentials_s{{simpl}}_{{clusters}}.csv",
                )
            },
        ),
        **branch(
            solar_rooftop_enable,
            {
                f"solar_rooftop_layout_{country}": "resources/"
                + RDIR
                + "solar_rooftop/solar_rooftop_layout_elec_s{simpl}_{clusters}_"
                + f"{country}.csv"
                for country in config["countries"]
            },
        ),
        network=rules.prepare_network.output.network,
        costs=rules.process_cost_data.output.costs.format(
            year="{planning_horizons}", scope="sec"
        ),
        nodal_energy_totals=branch(
            sector_enable["rail_transport"] or sector_enable["agriculture"],
            rules.prepare_heat_data.output.nodal_energy_totals,
        ),
        clustered_pop_layout=rules.build_clustered_population_layouts.output.clustered_pop_layout,
        industrial_demand=branch(
            sector_enable["industry"],
            rules.build_industry_demand.output.industrial_energy_demand_per_node,
        ),
        energy_totals=rules.prepare_energy_totals.output.energy_totals,
        airports=branch(
            sector_enable["aviation"],
            rules.prepare_airports.output.airports,
        ),
        ports=branch(sector_enable["shipping"], rules.prepare_ports.output.ports),
        biomass_transport_costs="data/temp_hard_coded/biomass_transport_costs.csv",
        shapes_path=rules.cluster_network.output.regions_onshore,
        pipelines=branch(
            config["sector"]["hydrogen"]["network"],
            branch(
                config["custom_data"]["gas_network"],
                "data/custom/pipelines.csv",
                "resources/"
                + SECDIR
                + "gas_networks/gas_network_elec_s{simpl}_{clusters}.csv",
            ),
        ),
    output:
        network=RESDIR
        + "prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}.nc",
    threads: 1
    resources:
        mem_mb=2000,
    benchmark:
        (
            RESDIR
            + "benchmarks/prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}"
        )
    script:
        scripts("prepare_sector_network.py")


rule build_ship_profile:
    params:
        snapshots=config["snapshots"],
        ship_opts=config["export"]["ship"],
        h2export=config["export"]["h2export"],
    output:
        ship_profile="resources/" + SECDIR + "ship_profile.csv",
    script:
        scripts("build_ship_profile.py")


rule add_export:
    params:
        gadm_layer_id=config["build_shape_options"]["gadm_layer_id"],
        alternative_clustering=config["clustering"]["alternative_clustering"],
        store=config["export"]["store"],
        store_capital_costs=config["export"]["store_capital_costs"],
        export_profile=config["export"]["export_profile"],
        export_endogenous=config["export"]["endogenous"],
        endogenous_price=config["export"]["endogenous_price"],
        snapshots=config["snapshots"],
        h2export=config["export"]["h2export"],
    input:
        export_ports=rules.prepare_ports.output.export_ports,
        costs=rules.process_cost_data.output.costs.format(
            year="{planning_horizons}", scope="sec"
        ),
        ship_profile=rules.build_ship_profile.output.ship_profile,
        network=rules.prepare_sector_network.output.network,
        shapes_path=rules.cluster_network.output.regions_onshore,
    output:
        network=RESDIR
        + "prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_export.nc",
    script:
        scripts("add_export.py")
