# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later


def retrieve_subregion(script_name):
    """
    Select whether scripts related to subregions should be retrieved.

    Offshore subregion shapes can be generated from subregion shapes, or provided as a separate custom file.
    """
    subregion_config = config.get("subregion", {"apply_on": []})

    if script_name not in subregion_config["apply_on"]:
        return {}

    method = subregion_config["method"]
    subregion_offshore = rules.build_shapes.output.offshore_shapes

    if method == "gadm":
        subregion_shapes = rules.build_shapes.output.subregion_shapes
    elif method == "custom":
        subregion_shapes = subregion_config["path_custom_shapes"]
        if subregion_config["path_custom_offshore"]:
            subregion_offshore = subregion_config["path_custom_offshore"]
    else:
        return {}

    return {
        "subregion_shapes": subregion_shapes,
        "subregion_offshore": subregion_offshore,
        "original_shapes": rules.build_shapes.output.country_shapes,
    }


rule build_shapes:
    params:
        build_shape_options=config["build_shape_options"],
        crs=config["crs"],
        countries=config["countries"],
        subregion=config["subregion"],
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
        extended_country_shape="resources/"
        + RDIR
        + "shapes/extended_country_shape.geojson",
        gadm_shapes="resources/" + RDIR + "shapes/gadm_shapes.geojson",
        subregion_shapes="resources/" + RDIR + "shapes/subregion_shapes.geojson",
        subregion_offshore="resources/" + RDIR + "shapes/subregion_offshore.geojson",
    log:
        "logs/" + RDIR + "build_shapes.log",
    benchmark:
        "benchmarks/" + RDIR + "build_shapes"
    threads: 1
    resources:
        mem_mb=3096,
    script:
        scripts("build_shapes.py")


rule clean_osm_data:
    params:
        crs=config["crs"],
        clean_osm_data_options=config["osm"]["clean_osm_data"],
    input:
        cables=rules.download_osm_data.output.cables,
        generators=rules.download_osm_data.output.generators,
        lines=rules.download_osm_data.output.lines,
        substations=rules.download_osm_data.output.substations,
        country_shapes=rules.build_shapes.output.country_shapes,
        offshore_shapes=rules.build_shapes.output.offshore_shapes,
        extended_country_shape=rules.build_shapes.output.extended_country_shape,
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
        scripts("clean_osm_data.py")


rule build_osm_network:
    params:
        build_osm_network=config.get("osm", {}).get("build_osm_network", {}),
        countries=config["countries"],
        crs=config["crs"],
    input:
        generators=rules.clean_osm_data.output.generators,
        lines=rules.clean_osm_data.output.lines,
        substations=rules.clean_osm_data.output.substations,
        country_shapes=rules.build_shapes.output.country_shapes,
    output:
        lines="resources/" + RDIR + "base_network/all_lines_build_network.csv",
        converters="resources/" + RDIR + "base_network/all_converters_build_network.csv",
        transformers="resources/"
        + RDIR
        + "base_network/all_transformers_build_network.csv",
        substations="resources/" + RDIR + "base_network/all_buses_build_network.csv",
        lines_geo="resources/" + RDIR + "base_network/all_lines_build_network.geojson",
        converters_geo="resources/"
        + RDIR
        + "base_network/all_converters_build_network.geojson",
        transformers_geo="resources/"
        + RDIR
        + "base_network/all_transformers_build_network.geojson",
        substations_geo="resources/"
        + RDIR
        + "base_network/all_buses_build_network.geojson",
    log:
        "logs/" + RDIR + "build_osm_network.log",
    benchmark:
        "benchmarks/" + RDIR + "build_osm_network"
    script:
        scripts("build_osm_network.py")


rule base_network:
    params:
        voltages=config["electricity"]["voltages"],
        transformers=config["transformers"],
        snapshots=config["snapshots"],
        links=config["links"],
        lines=config["lines"],
        hvdc_as_lines=config["electricity"]["hvdc_as_lines"],
        countries=config["countries"],
        base_network=config["base_network"],
    input:
        osm_buses=rules.build_osm_network.output.substations,
        osm_lines=rules.build_osm_network.output.lines,
        osm_converters=rules.build_osm_network.output.converters,
        osm_transformers=rules.build_osm_network.output.transformers,
        country_shapes=rules.build_shapes.output.country_shapes,
        offshore_shapes=rules.build_shapes.output.offshore_shapes,
        custom_line_types="data/custom_line_types.csv",
    output:
        network="networks/" + RDIR + "base.nc",
    log:
        "logs/" + RDIR + "base_network.log",
    benchmark:
        "benchmarks/" + RDIR + "base_network"
    threads: 1
    resources:
        mem_mb=500,
    script:
        scripts("base_network.py")


rule build_bus_regions:
    params:
        alternative_clustering=config["clustering"]["alternative_clustering"],
        crs=config["crs"],
        countries=config["countries"],
    input:
        **retrieve_subregion("cluster_network"),
        country_shapes=rules.build_shapes.output.country_shapes,
        offshore_shapes=rules.build_shapes.output.offshore_shapes,
        base_network=rules.base_network.output.network,
        #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson",
        #using this line instead of the following will test updated gadm shapes for MA.
        #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
        #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
        gadm_shapes=rules.build_shapes.output.gadm_shapes,
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
        scripts("build_bus_regions.py")


if config["enable"].get("build_cutout", False):

    rule build_cutout:
        params:
            snapshots=config["snapshots"],
            cutouts=config["atlite"]["cutouts"],
        input:
            check=terminate_if_cutout_exists,
            onshore_shapes=rules.build_shapes.output.country_shapes,
            offshore_shapes=rules.build_shapes.output.offshore_shapes,
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
            scripts("build_cutout.py")


if config["enable"].get("build_natura_raster", False):

    rule build_natura_raster:
        params:
            area_crs=config["crs"]["area_crs"],
            natura=config["natura"],
            disable_progress=not config["enable"]["progress_bar"],
        input:
            shapefiles_land="data/landcover",
            cutouts=expand(
                "cutouts/" + CDIR + "{cutout}.nc",
                cutout=[c["cutout"] for _, c in config["renewable"].items()],
            ),
            country_shapes=rules.build_shapes.output.country_shapes,
            offshore_shapes=rules.build_shapes.output.offshore_shapes,
        output:
            "resources/" + RDIR + "natura.tiff",
        log:
            "logs/" + RDIR + "build_natura_raster.log",
        benchmark:
            "benchmarks/" + RDIR + "build_natura_raster"
        script:
            scripts("build_natura_raster.py")


if not config["enable"].get("build_natura_raster", False):

    rule copy_defaultnatura_tiff:
        input:
            "data/natura/natura.tiff",
        output:
            "resources/" + RDIR + "natura.tiff",
        run:
            import shutil

            shutil.copyfile(input[0], output[0])


rule build_demand_profiles:
    params:
        snapshots=config["snapshots"],
        load_options=config["load_options"],
        countries=config["countries"],
    input:
        base_network=rules.base_network.output.network,
        regions=rules.build_bus_regions.output.regions_onshore,
        load=branch(
            config["load_options"].get("source", "gegis") in ["gegis", "ssp"],
            get_load_paths_gegis("data", config),
            "data/demand/forecasts_on_historical_period.parquet",
        ),
        #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson",
        #using this line instead of the following will test updated gadm shapes for MA.
        #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
        #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
        gadm_shapes=rules.build_shapes.output.gadm_shapes,
    output:
        demand_profiles="resources/" + RDIR + "demand_profiles.csv",
    log:
        "logs/" + RDIR + "build_demand_profiles.log",
    benchmark:
        "benchmarks/" + RDIR + "build_demand_profiles"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        scripts("build_demand_profiles.py")


def inputs_hydro(w):
    if w.technology == "hydro":
        HYDRO_PROFILES = {
            "hydro_capacities": "data/hydro_capacities.csv",
            "eia_hydro_generation": "data/eia_hydro_annual_generation.csv",
            "irena_stats": "data/IRENA_Statistics_Extract_2025H2.xlsx",
            "powerplants": rules.build_powerplants.output.powerplants,
            "hydrobasins": config["renewable"]["hydro"]["resource"]["hydrobasins"],
        }
        return HYDRO_PROFILES
    else:
        return {}


rule build_renewable_profiles:
    params:
        crs=config["crs"],
        renewable=config["renewable"],
        countries=config["countries"],
        alternative_clustering=config["clustering"]["alternative_clustering"],
    input:
        unpack(inputs_hydro),
        natura="resources/" + RDIR + "natura.tiff",
        copernicus="data/copernicus/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        gebco="data/gebco/GEBCO_2025_sub_ice.nc",
        country_shapes=rules.build_shapes.output.country_shapes,
        offshore_shapes=rules.build_shapes.output.offshore_shapes,
        regions=lambda w: (
            rules.build_bus_regions.output.regions_onshore
            if w.technology in ("onwind", "solar", "hydro", "csp")
            else rules.build_bus_regions.output.regions_offshore
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
        scripts("build_renewable_profiles.py")


rule build_powerplants:
    params:
        geo_crs=config["crs"]["geo_crs"],
        countries=config["countries"],
        gadm_layer_id=config["build_shape_options"]["gadm_layer_id"],
        alternative_clustering=config["clustering"]["alternative_clustering"],
        powerplants_filter=config["electricity"]["powerplants_filter"],
        custom_powerplants=config["electricity"]["custom_powerplants"],
    input:
        base_network=rules.base_network.output.network,
        pm_config="configs/powerplantmatching_config.yaml",
        custom_powerplants=branch(
            config["electricity"]["custom_powerplants"]["method"] is not False,
            config["electricity"]["custom_powerplants"]["filepaths"],
            [],
        ),
        osm_powerplants=rules.clean_osm_data.output.generators_csv,
        #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson",
        #using this line instead of the following will test updated gadm shapes for MA.
        #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
        #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
        gadm_shapes=rules.build_shapes.output.gadm_shapes,
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
        scripts("build_powerplants.py")


if config["enable"].get("retrieve_cost_data", True):
    cost_data = rules.retrieve_cost_data.output.costs
else:
    cost_data = "data/costs.csv"


rule process_cost_data:
    params:
        costs=config["costs"],
        max_hours=config["electricity"]["max_hours"],
        storage_techs=config["storage_techs"],
    input:
        network=rules.base_network.output.network,
        costs=cost_data,
    output:
        costs="resources/" + RDIR + "costs_{year}_{scope}.csv",
    log:
        "logs/" + RDIR + "build_cost_data_{year}_{scope}.log",
    benchmark:
        "benchmarks/" + RDIR + "build_cost_data_{year}_{scope}"
    threads: 1
    resources:
        mem_mb=4000,
    script:
        scripts("process_cost_data.py")


rule add_electricity:
    params:
        countries=config["countries"],
        output_currency=config["costs"]["output_currency"],
        fill_values=config["costs"]["fill_values"],
        conventional=config.get("conventional", {}),
        electricity=config["electricity"],
        alternative_clustering=config["clustering"]["alternative_clustering"],
        renewable=config["renewable"],
        length_factor=config["lines"]["length_factor"],
        existing_capacities=config["existing_capacities"],
        battery_techs=config["storage_techs"]["battery"],
    input:
        **{
            f"profile_{tech}": rules.build_renewable_profiles.output.profile.format(
                technology=tech
            )
            for tech in config["renewable"]
            if tech in config["electricity"]["renewable_carriers"]
        },
        **{
            f"conventional_{carrier}_{attr}": fn
            for carrier, d in config.get("conventional", {None: {}}).items()
            for attr, fn in d.items()
            if str(fn).startswith("data/")
        },
        base_network=rules.base_network.output.network,
        tech_costs=rules.process_cost_data.output.costs.format(
            year=config["costs"]["year"], scope="elec"
        ),
        powerplants=rules.build_powerplants.output.powerplants,
        #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson",
        #using this line instead of the following will test updated gadm shapes for MA.
        #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
        #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
        gadm_shapes=rules.build_shapes.output.gadm_shapes,
        hydro_capacities="data/hydro_capacities.csv",
        demand_profiles=rules.build_demand_profiles.output.demand_profiles,
        nuclear_p_max_pu="data/nuclear_p_max_pu.csv",
    output:
        network="networks/" + RDIR + "elec.nc",
    log:
        "logs/" + RDIR + "add_electricity.log",
    benchmark:
        "benchmarks/" + RDIR + "add_electricity"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        scripts("add_electricity.py")


rule simplify_network:
    params:
        aggregation_strategies=config["clustering"]["aggregation_strategies"],
        renewable=config["renewable"],
        crs=config["crs"],
        clustering=config["clustering"],
        countries=config["countries"],
        build_shape_options=config["build_shape_options"],
        electricity=config["electricity"],
        output_currency=config["costs"]["output_currency"],
        config_lines=config["lines"],
        config_links=config["links"],
        focus_weights=config.get("focus_weights", None),
    input:
        **retrieve_subregion("simplify_network"),
        network=rules.add_electricity.output.network,
        tech_costs=rules.process_cost_data.output.costs.format(
            year=config["costs"]["year"], scope="elec"
        ),
        regions_onshore=rules.build_bus_regions.output.regions_onshore,
        regions_offshore=rules.build_bus_regions.output.regions_offshore,
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
        scripts("simplify_network.py")


rule cluster_network:
    params:
        aggregation_strategies=config["clustering"]["aggregation_strategies"],
        build_shape_options=config["build_shape_options"],
        electricity=config["electricity"],
        length_factor=config["lines"]["length_factor"],
        renewable=config["renewable"],
        crs=config["crs"],
        countries=config["countries"],
        clustering=config["clustering"],
        focus_weights=config.get("focus_weights", None),
        custom_busmap=config["enable"].get("custom_busmap", False),
    input:
        **retrieve_subregion("cluster_network"),
        network=rules.simplify_network.output.network,
        country_shapes=rules.build_shapes.output.country_shapes,
        regions_onshore=rules.simplify_network.output.regions_onshore,
        regions_offshore=rules.simplify_network.output.regions_offshore,
        #using this line instead of the following will test updated gadm shapes for MA.
        #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
        #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
        gadm_shapes=rules.build_shapes.output.gadm_shapes,
        # busmap=ancient('resources/" + RDIR + "bus_regions/busmap_elec_s{simpl}.csv'),
        custom_busmap=(
            "data/custom_busmap_elec_s{simpl}_{clusters}.csv"
            if config["enable"].get("custom_busmap", False)
            else []
        ),
        tech_costs=rules.process_cost_data.output.costs.format(
            year=config["costs"]["year"], scope="elec"
        ),
    output:
        network=branch(
            config["augmented_line_connection"].get("add_to_snakefile", False) == True,
            "networks/" + RDIR + "elec_s{simpl}_{clusters}_pre_augmentation.nc",
            "networks/" + RDIR + "elec_s{simpl}_{clusters}.nc",
        ),
        regions_onshore="resources/"
        + RDIR
        + "bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
        regions_offshore="resources/"
        + RDIR
        + "bus_regions/regions_offshore_elec_s{simpl}_{clusters}.geojson",
        busmap="resources/" + RDIR + "bus_regions/busmap_elec_s{simpl}_{clusters}.csv",
        linemap="resources/" + RDIR + "bus_regions/linemap_elec_s{simpl}_{clusters}.csv",
    log:
        "logs/" + RDIR + "cluster_network/elec_s{simpl}_{clusters}.log",
    benchmark:
        "benchmarks/" + RDIR + "cluster_network/elec_s{simpl}_{clusters}"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        scripts("cluster_network.py")


if config["augmented_line_connection"].get("add_to_snakefile") == True:

    rule augmented_line_connections:
        params:
            lines=config["lines"],
            augmented_line_connection=config["augmented_line_connection"],
            hvdc_as_lines=config["electricity"]["hvdc_as_lines"],
            electricity=config["electricity"],
        input:
            tech_costs=rules.process_cost_data.output.costs.format(
                year=config["costs"]["year"], scope="elec"
            ),
            network=rules.cluster_network.output.network,
            regions_onshore=rules.cluster_network.output.regions_onshore,
            regions_offshore=rules.cluster_network.output.regions_offshore,
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
            scripts("augmented_line_connections.py")


if config["augmented_line_connection"].get("add_to_snakefile") == True:
    clustered_network = rules.augmented_line_connections.output.network
else:
    clustered_network = rules.cluster_network.output.network


rule add_extra_components:
    params:
        storage_techs=config["storage_techs"],
        transmission_efficiency=config["sector"]["transmission_efficiency"],
        electricity=config["electricity"],
        csp_model=config["renewable"]["csp"]["csp_model"],
    input:
        network=clustered_network,
        tech_costs=rules.process_cost_data.output.costs.format(
            year=config["costs"]["year"], scope="elec"
        ),
    output:
        network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec.nc",
    log:
        "logs/" + RDIR + "add_extra_components/elec_s{simpl}_{clusters}.log",
    benchmark:
        "benchmarks/" + RDIR + "add_extra_components/elec_s{simpl}_{clusters}_ec"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        scripts("add_extra_components.py")


if config["co2"]["automatic_emission"]["enable"]:

    rule build_co2_emissions:
        input:
            edgar=rules.retrieve_emissions.output.edgar_xlsx,
        output:
            emissions="resources/" + RDIR + "co2_emissions_elec_and_heat.csv",
        log:
            "logs/" + RDIR + "build_co2_emissions.log",
        benchmark:
            "benchmarks/" + RDIR + "build_co2_emissions"
        resources:
            mem_mb=2000,
        script:
            scripts("build_co2_emissions.py")


if config["co2"]["automatic_emission"]["enable"]:
    emissions_input = {"emissions": rules.build_co2_emissions.output.emissions}
else:
    emissions_input = {}


rule prepare_network:
    params:
        links=config["links"],
        lines=config["lines"],
        s_max_pu=config["lines"]["s_max_pu"],
        electricity=config["electricity"],
        co2=config["co2"],
    input:
        **emissions_input,
        network=rules.add_extra_components.output.network,
        tech_costs=rules.process_cost_data.output.costs.format(
            year=config["costs"]["year"], scope="elec"
        ),
    output:
        network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
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
        scripts("prepare_network.py")


if config["monte_carlo"]["options"].get("add_to_snakefile", False) == True:

    rule monte_carlo:
        params:
            monte_carlo=config["monte_carlo"],
        input:
            network=rules.prepare_network.output.network,
        output:
            network="networks/"
            + RDIR
            + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}.nc",
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
            scripts("monte_carlo.py")
