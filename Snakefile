import sys


PYPSAEARTH_FOLDER = "pypsa-earth"
sys.path.append(PYPSAEARTH_FOLDER + "/scripts")

from os.path import exists
from pathlib import Path
from shutil import copyfile, move
from scripts.helpers import get_last_commit_message

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from _helpers import create_country_list

HTTP = HTTPRemoteProvider()

if not exists("config.yaml"):
    copyfile("config.default.yaml", "config.yaml")


configfile: "config.pypsa-earth.yaml"
configfile: PYPSAEARTH_FOLDER + "/configs/bundle_config.yaml"
configfile: PYPSAEARTH_FOLDER + "/configs/powerplantmatching_config.yaml"
configfile: "config.yaml"


# convert country list according to the desired region
config["countries"] = create_country_list(config["countries"])
config.update({"git_commit": get_last_commit_message(".")})
config.update({"submodule_commit": get_last_commit_message(PYPSAEARTH_FOLDER)})
config["ROOT_PATH"] = os.getcwd()

run = config.get("run", {})
RDIR = run["name"] + "/" if run.get("name") else ""
CDIR = RDIR if not run.get("shared_cutouts") else ""
SDIR = config["summary_dir"] + (run["name"] + "/" if run.get("name") else "")
RESDIR = config["results_dir"] + (run["name"] + "/" if run.get("name") else "")
COSTDIR = config["costs_dir"]

# RDIR_PE = run["name_subworkflow"] + "/" if run.get("name_subworkflow") else ""
# CDIR_PE = RDIR_PE if not run.get("shared_cutouts") else ""


CUTOUTS_PATH = (
    "cutouts/"
    + CDIR
    + ("cutout-2013-era5-tutorial.nc" if config["tutorial"] else "cutout-2013-era5.nc")
)


wildcard_constraints:
    ll="[a-z0-9\.]+",
    simpl="[a-zA-Z0-9]*|all",
    clusters="[0-9]+m?|all",
    opts="[-+a-zA-Z0-9\.\s]*",
    sopts="[-+a-zA-Z0-9\.\s]*",
    discountrate="[-+a-zA-Z0-9\.\s]*",
    demand="[-+a-zA-Z0-9\.\s]*",
    h2export="[0-9]+m?|all",
    planning_horizons="20[2-9][0-9]|2100",


if not config.get("disable_subworkflow", False):

    module pypsaearth:
        snakefile:
            PYPSAEARTH_FOLDER + "/Snakefile"
        config:
            config

    use rule * from pypsaearth


data_dir = Path(PYPSAEARTH_FOLDER) / "data"


rule get_data:
    output:
        [
            str(Path("data") / p.relative_to(data_dir))
            for p in data_dir.rglob("*")
            if p.is_file()
        ],
    shell:
        """
        mkdir -p data
        cp -nR {data_dir}/. data/
        """


if config["enable"].get("retrieve_cost_data", True):

    rule retrieve_cost_data_flexible:
        input:
            HTTP.remote(
                f"raw.githubusercontent.com/PyPSA/technology-data/{config['costs']['version']}/outputs/costs"
                + "_{planning_horizons}.csv",
                keep_local=True,
            ),
        output:
            costs=COSTDIR + "costs_{planning_horizons}.csv",
        resources:
            mem_mb=5000,
        run:
            move(input[0], output[0])


rule prepare_sector_networks:
    input:
        expand(
            RESDIR
            + "prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}.nc",
            **config["scenario"],
            **config["costs"],
        ),


rule override_res_all_nets:
    input:
        expand(
            RESDIR
            + "prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_presec.nc",
            **config["scenario"],
            **config["costs"],
            **config["export"],
        ),


rule solve_sector_networks:
    input:
        expand(
            RESDIR
            + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
            **config["scenario"],
            **config["costs"],
            **config["export"],
        ),


rule prepare_ports:
    output:
        ports="data/ports.csv",  # TODO move from data to resources
    script:
        "scripts/prepare_ports.py"


rule prepare_airports:
    params:
        airport_sizing_factor=config["sector"]["airport_sizing_factor"],
    output:
        ports="data/airports.csv",  # TODO move from data to resources
    script:
        "scripts/prepare_airports.py"


rule prepare_urban_percent:
    output:
        urban_percent="data/urban_percent.csv",  # TODO move from data to resources
    script:
        "scripts/prepare_urban_percent.py"


rule prepare_transport_data_input:
    output:
        transport_data_input="resources/transport_data.csv",
    script:
        "scripts/prepare_transport_data_input.py"


if not config["custom_data"]["gas_network"]:

    rule prepare_gas_network:
        params:
            gas_config=config["sector"]["gas"],
            alternative_clustering=config["clustering_options"][
                "alternative_clustering"
            ],
            countries_list=config["countries"],
            layer_id=config["build_shape_options"]["gadm_layer_id"],
            update=config["build_shape_options"]["update_file"],
            out_logging=config["build_shape_options"]["out_logging"],
            year=config["build_shape_options"]["year"],
            nprocesses=config["build_shape_options"]["nprocesses"],
            contended_flag=config["build_shape_options"]["contended_flag"],
            geo_crs=config["crs"]["geo_crs"],
            custom_gas_network=config["custom_data"]["gas_network"],
        input:
            regions_onshore="resources/"
            + RDIR
            + "bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
        output:
            clustered_gas_network="resources/gas_networks/gas_network_elec_s{simpl}_{clusters}.csv",
            # TODO: Should be a own snakemake rule
            # gas_network_fig_1="resources/gas_networks/existing_gas_pipelines_{simpl}_{clusters}.png",
            # gas_network_fig_2="resources/gas_networks/clustered_gas_pipelines_{simpl}_{clusters}.png",
        script:
            "scripts/prepare_gas_network.py"


rule prepare_sector_network:
    input:
        network=RESDIR
        + "prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_presec.nc",
        costs=COSTDIR + "costs_{planning_horizons}.csv",
        h2_cavern="data/hydrogen_salt_cavern_potentials.csv",
        nodal_energy_totals="resources/demand/heat/nodal_energy_heat_totals_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        transport="resources/demand/transport_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        avail_profile="resources/pattern_profiles/avail_profile_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        dsm_profile="resources/pattern_profiles/dsm_profile_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        nodal_transport_data="resources/demand/nodal_transport_data_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        overrides="data/override_component_attrs",
        clustered_pop_layout="resources/population_shares/pop_layout_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
        industrial_demand="resources/demand/industrial_energy_demand_per_node_elec_s{simpl}_{clusters}_{planning_horizons}_{demand}.csv",
        energy_totals="data/energy_totals_{demand}_{planning_horizons}.csv",
        airports="data/airports.csv",
        ports="data/ports.csv",
        heat_demand="resources/demand/heat/heat_demand_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        ashp_cop="resources/demand/heat/ashp_cop_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        gshp_cop="resources/demand/heat/gshp_cop_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        solar_thermal="resources/demand/heat/solar_thermal_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        district_heat_share="resources/demand/heat/district_heat_share_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        biomass_transport_costs="data/temp_hard_coded/biomass_transport_costs.csv",
        shapes_path="resources/"
        + RDIR
        + "bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
        pipelines=(
            "data_custom/pipelines.csv"
            if config["custom_data"]["gas_network"]
            else "resources/gas_networks/gas_network_elec_s{simpl}_{clusters}.csv"
        ),
    output:
        RESDIR
        + "prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}.nc",
    threads: 1
    resources:
        mem_mb=2000,
    benchmark:
        (
            RESDIR
            + "/benchmarks/prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}"
        )
    script:
        "scripts/prepare_sector_network.py"


rule build_ship_profile:
    params:
        snapshots=config["snapshots"],
        ship_opts=config["export"]["ship"],
    output:
        ship_profile="resources/ship_profile_{h2export}TWh.csv",
    script:
        "scripts/build_ship_profile.py"


rule add_export:
    params:
        gadm_level=config["sector"]["gadm_level"],
        alternative_clustering=config["clustering_options"]["alternative_clustering"],
        store=config["export"]["store"],
        store_capital_costs=config["export"]["store_capital_costs"],
        export_profile=config["export"]["export_profile"],
        snapshots=config["snapshots"],
        USD_to_EUR=config["costs"]["USD2013_to_EUR2013"],
        lifetime=config["costs"]["lifetime"],
    input:
        overrides="data/override_component_attrs",
        export_ports="data/export_ports.csv",
        costs=COSTDIR + "costs_{planning_horizons}.csv",
        ship_profile="resources/ship_profile_{h2export}TWh.csv",
        network=RESDIR
        + "prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}.nc",
        shapes_path="resources/"
        + RDIR
        + "bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
    output:
        RESDIR
        + "prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
    script:
        "scripts/add_export.py"


rule override_respot:
    params:
        run=run["name"],
        custom_data=config["custom_data"],
        countries=config["countries"],
    input:
        **{
            f"custom_res_pot_{tech}_{planning_horizons}_{discountrate}": f"resources/custom_renewables/{tech}_{planning_horizons}_{discountrate}_potential.csv"
            for tech in config["custom_data"]["renewables"]
            for discountrate in config["costs"]["discountrate"]
            for planning_horizons in config["scenario"]["planning_horizons"]
        },
        **{
            f"custom_res_ins_{tech}_{planning_horizons}_{discountrate}": f"resources/custom_renewables/{tech}_{planning_horizons}_{discountrate}_installable.csv"
            for tech in config["custom_data"]["renewables"]
            for discountrate in config["costs"]["discountrate"]
            for planning_horizons in config["scenario"]["planning_horizons"]
        },
        overrides="data/override_component_attrs",
        network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        energy_totals="data/energy_totals_{demand}_{planning_horizons}.csv",
    output:
        RESDIR
        + "prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_presec.nc",
    script:
        "scripts/override_respot.py"


rule prepare_transport_data:
    input:
        network="networks/" + RDIR + "elec_s{simpl}_{clusters}.nc",
        energy_totals_name="data/energy_totals_{demand}_{planning_horizons}.csv",
        traffic_data_KFZ="data/emobility/KFZ__count",
        traffic_data_Pkw="data/emobility/Pkw__count",
        transport_name="resources/transport_data.csv",
        clustered_pop_layout="resources/population_shares/pop_layout_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
        temp_air_total="resources/temperatures/temp_air_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
    output:
        # nodal_energy_totals="resources/nodal_energy_totals_s{simpl}_{clusters}.csv",
        transport="resources/demand/transport_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        avail_profile="resources/pattern_profiles/avail_profile_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        dsm_profile="resources/pattern_profiles/dsm_profile_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        nodal_transport_data="resources/demand/nodal_transport_data_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
    script:
        "scripts/prepare_transport_data.py"


rule build_cop_profiles:
    params:
        heat_pump_sink_T=config["sector"]["heat_pump_sink_T"],
    input:
        temp_soil_total="resources/temperatures/temp_soil_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        temp_soil_rural="resources/temperatures/temp_soil_rural_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        temp_soil_urban="resources/temperatures/temp_soil_urban_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        temp_air_total="resources/temperatures/temp_air_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        temp_air_rural="resources/temperatures/temp_air_rural_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        temp_air_urban="resources/temperatures/temp_air_urban_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
    output:
        cop_soil_total="resources/cops/cop_soil_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        cop_soil_rural="resources/cops/cop_soil_rural_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        cop_soil_urban="resources/cops/cop_soil_urban_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        cop_air_total="resources/cops/cop_air_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        cop_air_rural="resources/cops/cop_air_rural_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        cop_air_urban="resources/cops/cop_air_urban_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        "benchmarks/build_cop_profiles/s{simpl}_{clusters}_{planning_horizons}"
    script:
        "scripts/build_cop_profiles.py"


rule prepare_heat_data:
    input:
        network="networks/" + RDIR + "elec_s{simpl}_{clusters}.nc",
        energy_totals_name="data/energy_totals_{demand}_{planning_horizons}.csv",
        clustered_pop_layout="resources/population_shares/pop_layout_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
        temp_air_total="resources/temperatures/temp_air_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        cop_soil_total="resources/cops/cop_soil_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        cop_air_total="resources/cops/cop_air_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        solar_thermal_total="resources/demand/heat/solar_thermal_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        heat_demand_total="resources/demand/heat/heat_demand_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        heat_profile="data/heat_load_profile_BDEW.csv",
    output:
        nodal_energy_totals="resources/demand/heat/nodal_energy_heat_totals_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        heat_demand="resources/demand/heat/heat_demand_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        ashp_cop="resources/demand/heat/ashp_cop_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        gshp_cop="resources/demand/heat/gshp_cop_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        solar_thermal="resources/demand/heat/solar_thermal_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        district_heat_share="resources/demand/heat/district_heat_share_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
    script:
        "scripts/prepare_heat_data.py"


rule build_base_energy_totals:
    params:
        space_heat_share=config["sector"]["space_heat_share"],
        update_data=config["demand_data"]["update_data"],
        base_year=config["demand_data"]["base_year"],
        countries=config["countries"],
    input:
        unsd_paths="data/demand/unsd/paths/Energy_Statistics_Database.xlsx",
    output:
        energy_totals_base="data/energy_totals_base.csv",
    script:
        "scripts/build_base_energy_totals.py"


rule prepare_energy_totals:
    params:
        countries=config["countries"],
        base_year=config["demand_data"]["base_year"],
        sector_options=config["sector"],
    input:
        unsd_paths="data/energy_totals_base.csv",
        efficiency_gains_cagr="data/demand/efficiency_gains_cagr.csv",
        growth_factors_cagr="data/demand/growth_factors_cagr.csv",
        district_heating="data/demand/district_heating.csv",
        fuel_shares="data/demand/fuel_shares.csv",
    output:
        energy_totals="data/energy_totals_{demand}_{planning_horizons}.csv",
    script:
        "scripts/prepare_energy_totals.py"


rule build_solar_thermal_profiles:
    params:
        solar_thermal_config=config["solar_thermal"],
        snapshots=config["snapshots"],
    input:
        pop_layout_total="resources/population_shares/pop_layout_total_{planning_horizons}.nc",
        pop_layout_urban="resources/population_shares/pop_layout_urban_{planning_horizons}.nc",
        pop_layout_rural="resources/population_shares/pop_layout_rural_{planning_horizons}.nc",
        regions_onshore="resources/"
        + RDIR
        + "bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
        cutout=CUTOUTS_PATH,
    output:
        solar_thermal_total="resources/demand/heat/solar_thermal_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        solar_thermal_urban="resources/demand/heat/solar_thermal_urban_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        solar_thermal_rural="resources/demand/heat/solar_thermal_rural_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        "benchmarks/build_solar_thermal_profiles/s{simpl}_{clusters}_{planning_horizons}"
    script:
        "scripts/build_solar_thermal_profiles.py"


rule build_population_layouts:
    params:
        planning_horizons=config["scenario"]["planning_horizons"][0],
    input:
        nuts3_shapes="resources/" + RDIR + "shapes/gadm_shapes.geojson",
        urban_percent="data/urban_percent.csv",
        cutout=CUTOUTS_PATH,
    output:
        pop_layout_total="resources/population_shares/pop_layout_total_{planning_horizons}.nc",
        pop_layout_urban="resources/population_shares/pop_layout_urban_{planning_horizons}.nc",
        pop_layout_rural="resources/population_shares/pop_layout_rural_{planning_horizons}.nc",
        gdp_layout="resources/gdp_shares/gdp_layout_{planning_horizons}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        "benchmarks/build_population_layouts_{planning_horizons}"
    threads: 8
    script:
        "scripts/build_population_layouts.py"


rule move_hardcoded_files_temp:
    input:
        "data/temp_hard_coded/energy_totals.csv",
    output:
        "resources/energy_totals.csv",
    shell:
        "cp -a data/temp_hard_coded/. resources"


rule build_clustered_population_layouts:
    input:
        pop_layout_total="resources/population_shares/pop_layout_total_{planning_horizons}.nc",
        pop_layout_urban="resources/population_shares/pop_layout_urban_{planning_horizons}.nc",
        pop_layout_rural="resources/population_shares/pop_layout_rural_{planning_horizons}.nc",
        gdp_layout="resources/gdp_shares/gdp_layout_{planning_horizons}.nc",
        regions_onshore="resources/"
        + RDIR
        + "bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
        cutout=CUTOUTS_PATH,
    output:
        clustered_pop_layout="resources/population_shares/pop_layout_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
        clustered_gdp_layout="resources/gdp_shares/gdp_layout_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
    resources:
        mem_mb=10000,
    benchmark:
        "benchmarks/build_clustered_population_layouts/s{simpl}_{clusters}_{planning_horizons}"
    script:
        "scripts/build_clustered_population_layouts.py"


rule build_heat_demand:
    params:
        snapshots=config["snapshots"],
    input:
        pop_layout_total="resources/population_shares/pop_layout_total_{planning_horizons}.nc",
        pop_layout_urban="resources/population_shares/pop_layout_urban_{planning_horizons}.nc",
        pop_layout_rural="resources/population_shares/pop_layout_rural_{planning_horizons}.nc",
        regions_onshore="resources/"
        + RDIR
        + "bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
        cutout=CUTOUTS_PATH,
    output:
        heat_demand_urban="resources/demand/heat/heat_demand_urban_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        heat_demand_rural="resources/demand/heat/heat_demand_rural_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        heat_demand_total="resources/demand/heat/heat_demand_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        "benchmarks/build_heat_demand/s{simpl}_{clusters}_{planning_horizons}"
    script:
        "scripts/build_heat_demand.py"


rule build_temperature_profiles:
    params:
        snapshots=config["snapshots"],
    input:
        pop_layout_total="resources/population_shares/pop_layout_total_{planning_horizons}.nc",
        pop_layout_urban="resources/population_shares/pop_layout_urban_{planning_horizons}.nc",
        pop_layout_rural="resources/population_shares/pop_layout_rural_{planning_horizons}.nc",
        regions_onshore="resources/"
        + RDIR
        + "bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
        cutout=CUTOUTS_PATH,
    output:
        temp_soil_total="resources/temperatures/temp_soil_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        temp_soil_rural="resources/temperatures/temp_soil_rural_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        temp_soil_urban="resources/temperatures/temp_soil_urban_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        temp_air_total="resources/temperatures/temp_air_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        temp_air_rural="resources/temperatures/temp_air_rural_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        temp_air_urban="resources/temperatures/temp_air_urban_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        "benchmarks/build_temperature_profiles/s{simpl}_{clusters}_{planning_horizons}"
    script:
        "scripts/build_temperature_profiles.py"


rule copy_config:
    params:
        summary_dir=config["summary_dir"],
        run=run,
    output:
        folder=directory(SDIR + "/configs"),
        config=SDIR + "/configs/config.yaml",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        SDIR + "/benchmarks/copy_config"
    script:
        "scripts/copy_config.py"


if config["foresight"] == "overnight":

    rule solve_sector_network:
        input:
            overrides="data/override_component_attrs",
            # network=RESDIR
            # + "prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}.nc",
            network=RESDIR
            + "prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
            costs=COSTDIR + "costs_{planning_horizons}.csv",
            configs=SDIR + "/configs/config.yaml",  # included to trigger copy_config rule
        output:
            RESDIR
            + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
        shadow:
            "shallow"
        log:
            solver=RESDIR
            + "/logs/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export_solver.log",
            python=RESDIR
            + "/logs/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export_python.log",
            memory=RESDIR
            + "/logs/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export_memory.log",
        threads: 25
        resources:
            mem_mb=config["solving"]["mem"],
        benchmark:
            (
                RESDIR
                + "/benchmarks/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export"
            )
        script:
            "scripts/solve_network.py"


rule make_sector_summary:
    params:
        planning_horizons=config["scenario"]["planning_horizons"],
        results_dir=config["results_dir"],
        summary_dir=config["summary_dir"],
        run=run["name"],
        scenario_config=config["scenario"],
        costs_config=config["costs"],
        h2export_qty=config["export"]["h2export"],
        foresight=config["foresight"],
    input:
        overrides="data/override_component_attrs",
        networks=expand(
            RESDIR
            + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
            **config["scenario"],
            **config["costs"],
            **config["export"],
        ),
        costs=COSTDIR + "costs_{planning_horizons}.csv",
        plots=expand(
            RESDIR
            + "/maps/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}-costs-all_{planning_horizons}_{discountrate}_{demand}_{h2export}export.pdf",
            **config["scenario"],
            **config["costs"],
            **config["export"],
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


rule plot_sector_network:
    input:
        overrides="data/override_component_attrs",
        network=RESDIR
        + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
    output:
        map=RESDIR
        + "/maps/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}-costs-all_{planning_horizons}_{discountrate}_{demand}_{h2export}export.pdf",
    threads: 2
    resources:
        mem_mb=10000,
    benchmark:
        (
            RESDIR
            + "/benchmarks/plot_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export"
        )
    script:
        "scripts/plot_network.py"


rule plot_sector_summary:
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


rule build_industrial_database:
    output:
        industrial_database="data/industrial_database.csv",
    script:
        "scripts/build_industrial_database.py"


rule prepare_db:
    params:
        tech_colors=config["plotting"]["tech_colors"],
    input:
        network=RESDIR
        + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
    output:
        db=RESDIR
        + "/summaries/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}-costs-all_{planning_horizons}_{discountrate}_{demand}_{h2export}export.csv",
    threads: 2
    resources:
        mem_mb=10000,
    benchmark:
        (
            RESDIR
            + "/benchmarks/prepare_db/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export"
        )
    script:
        "scripts/prepare_db.py"


rule build_industrial_distribution_key:  #default data
    params:
        countries=config["countries"],
        gadm_level=config["sector"]["gadm_level"],
        alternative_clustering=config["clustering_options"]["alternative_clustering"],
        industry_database=config["custom_data"]["industry_database"],
    input:
        regions_onshore="resources/"
        + RDIR
        + "bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
        clustered_pop_layout="resources/population_shares/pop_layout_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
        clustered_gdp_layout="resources/gdp_shares/gdp_layout_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
        industrial_database="data/industrial_database.csv",
        shapes_path="resources/"
        + RDIR
        + "bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
    output:
        industrial_distribution_key="resources/demand/industrial_distribution_key_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/build_industrial_distribution_key_elec_s{simpl}_{clusters}_{planning_horizons}"
    script:
        "scripts/build_industrial_distribution_key.py"


rule build_base_industry_totals:  #default data
    params:
        base_year=config["demand_data"]["base_year"],
        countries=config["countries"],
        other_industries=config["demand_data"]["other_industries"],
    input:
        #industrial_production_per_country="data/industrial_production_per_country.csv",
        #unsd_path="data/demand/unsd/data/",
        energy_totals_base="data/energy_totals_base.csv",
    output:
        base_industry_totals="resources/demand/base_industry_totals_{planning_horizons}_{demand}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/build_base_industry_totals_{planning_horizons}_{demand}"
    script:
        "scripts/build_base_industry_totals.py"


rule build_industry_demand:  #default data
    params:
        countries=config["countries"],
        industry_demand=config["custom_data"]["industry_demand"],
        base_year=config["demand_data"]["base_year"],
        industry_util_factor=config["sector"]["industry_util_factor"],
        aluminium_year=config["demand_data"]["aluminium_year"],
    input:
        industrial_distribution_key="resources/demand/industrial_distribution_key_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
        #industrial_production_per_country_tomorrow="resources/demand/industrial_production_per_country_tomorrow_{planning_horizons}_{demand}.csv",
        #industrial_production_per_country="data/industrial_production_per_country.csv",
        base_industry_totals="resources/demand/base_industry_totals_{planning_horizons}_{demand}.csv",
        industrial_database="data/industrial_database.csv",
        costs=COSTDIR + "costs_{planning_horizons}.csv",
        industry_growth_cagr="data/demand/industry_growth_cagr.csv",
    output:
        industrial_energy_demand_per_node="resources/demand/industrial_energy_demand_per_node_elec_s{simpl}_{clusters}_{planning_horizons}_{demand}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/industrial_energy_demand_per_node_elec_s{simpl}_{clusters}_{planning_horizons}_{demand}.csv"
    script:
        "scripts/build_industry_demand.py"


rule build_existing_heating_distribution:
    params:
        baseyear=config["scenario"]["planning_horizons"][0],
        sector=config["sector"],
        existing_capacities=config["existing_capacities"],
    input:
        existing_heating="data/existing_infrastructure/existing_heating_raw.csv",
        clustered_pop_layout="resources/population_shares/pop_layout_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
        clustered_pop_energy_layout="resources/demand/heat/nodal_energy_heat_totals_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",  #"resources/population_shares/pop_weighted_energy_totals_s{simpl}_{clusters}.csv",
        district_heat_share="resources/demand/heat/district_heat_share_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
    output:
        existing_heating_distribution="resources/heating/existing_heating_distribution_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
    threads: 1
    resources:
        mem_mb=2000,
    log:
        RESDIR
        + "/logs/build_existing_heating_distribution_{demand}_s{simpl}_{clusters}_{planning_horizons}.log",
    benchmark:
        RESDIR
        +"/benchmarks/build_existing_heating_distribution/{demand}_s{simpl}_{clusters}_{planning_horizons}"
    script:
        "scripts/build_existing_heating_distribution.py"


if config["foresight"] == "myopic":

    rule add_existing_baseyear:
        params:
            baseyear=config["scenario"]["planning_horizons"][0],
            sector=config["sector"],
            existing_capacities=config["existing_capacities"],
            costs=config["costs"],
        input:
            network=RESDIR
            + "/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
            powerplants="resources/" + RDIR + "powerplants.csv",
            busmap_s="resources/" + RDIR + "bus_regions/busmap_elec_s{simpl}.csv",
            busmap=pypsaearth(
                "resources/" + RDIR + "bus_regions/busmap_elec_s{simpl}_{clusters}.csv"
            ),
            clustered_pop_layout="resources/population_shares/pop_layout_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
            costs=CDIR
            + "costs_{}.csv".format(config["scenario"]["planning_horizons"][0]),
            cop_soil_total="resources/cops/cop_soil_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
            cop_air_total="resources/cops/cop_air_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
            existing_heating_distribution="resources/heating/existing_heating_distribution_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        output:
            RESDIR
            + "/prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
        wildcard_constraints:
            # TODO: The first planning_horizon needs to be aligned across scenarios
            # snakemake does not support passing functions to wildcard_constraints
            # reference: https://github.com/snakemake/snakemake/issues/2703
            planning_horizons=config["scenario"]["planning_horizons"][0],  #only applies to baseyear
        threads: 1
        resources:
            mem_mb=2000,
        log:
            RESDIR
            + "/logs/add_existing_baseyear_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.log",
        benchmark:
            RESDIR
            +"/benchmarks/add_existing_baseyear/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export"
        script:
            "scripts/add_existing_baseyear.py"

    def input_profile_tech_brownfield(w):
        return {
            f"profile_{tech}": f"resources/"
            + RDIR
            + "renewable_profiles/profile_{tech}.nc"
            for tech in config["electricity"]["renewable_carriers"]
            if tech != "hydro"
        }

    def solved_previous_horizon(w):
        planning_horizons = config["scenario"]["planning_horizons"]
        i = planning_horizons.index(int(w.planning_horizons))
        planning_horizon_p = str(planning_horizons[i - 1])

        return (
            RDIR
            + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_"
            + planning_horizon_p
            + "_{discountrate}_{demand}_{h2export}export.nc"
        )

    rule add_brownfield:
        params:
            H2_retrofit=config["sector"]["hydrogen"],
            H2_retrofit_capacity_per_CH4=config["sector"]["hydrogen"][
                "H2_retrofit_capacity_per_CH4"
            ],
            threshold_capacity=config["existing_capacities"]["threshold_capacity"],
            snapshots=config["snapshots"],
            # drop_leap_day=config["enable"]["drop_leap_day"],
            carriers=config["electricity"]["renewable_carriers"],
        input:
            # unpack(input_profile_tech_brownfield),
            simplify_busmap="resources/" + RDIR + "bus_regions/busmap_elec_s{simpl}.csv",
            cluster_busmap="resources/"
            + RDIR
            + "bus_regions/busmap_elec_s{simpl}_{clusters}.csv",
            network=RESDIR
            + "/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
            network_p=solved_previous_horizon,  #solved network at previous time step
            costs=CDIR + "costs_{planning_horizons}.csv",
            cop_soil_total="resources/cops/cop_soil_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
            cop_air_total="resources/cops/cop_air_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc",
        output:
            RESDIR
            + "/prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
        threads: 4
        resources:
            mem_mb=10000,
        log:
            RESDIR
            + "/logs/add_brownfield_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.log",
        benchmark:
            (
                RESDIR
                + "/benchmarks/add_brownfield/elec_s{simpl}_ec_{clusters}_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export"
            )
        script:
            "./scripts/add_brownfield.py"

    ruleorder: add_existing_baseyear > add_brownfield

    rule solve_network_myopic:
        params:
            solving=config["solving"],
            foresight=config["foresight"],
            planning_horizons=config["scenario"]["planning_horizons"],
            co2_sequestration_potential=config["scenario"].get(
                "co2_sequestration_potential", 200
            ),
        input:
            overrides="data/override_component_attrs",
            network=RESDIR
            + "/prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
            costs=CDIR + "costs_{planning_horizons}.csv",
            configs=SDIR + "/configs/config.yaml",  # included to trigger copy_config rule
        output:
            network=RESDIR
            + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
            # config=RESDIR
            # + "/configs/config.elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.yaml",
        shadow:
            "shallow"
        log:
            solver=RESDIR
            + "/logs/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export_solver.log",
            python=RESDIR
            + "/logs/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export_python.log",
            memory=RESDIR
            + "/logs/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export_memory.log",
        threads: 25
        resources:
            mem_mb=config["solving"]["mem"],
        benchmark:
            (
                RESDIR
                + "/benchmarks/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export"
            )
        script:
            "./scripts/solve_network.py"

    rule solve_all_networks_myopic:
        input:
            networks=expand(
                RESDIR
                + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
                **config["scenario"],
                **config["costs"],
                **config["export"],
            ),
