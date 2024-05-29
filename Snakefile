import sys

sys.path.append("pypsa-earth/scripts")

from os.path import exists
from shutil import copyfile, move
from scripts.helpers import get_last_commit_message

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from _helpers import create_country_list

HTTP = HTTPRemoteProvider()

if not exists("config.yaml"):
    copyfile("config.default.yaml", "config.yaml")


configfile: "config.pypsa-earth.yaml"
configfile: "config.yaml"


PYPSAEARTH_FOLDER = "pypsa-earth"

# convert country list according to the desired region
config["countries"] = create_country_list(config["countries"])


SDIR = config["summary_dir"] + config["run"]
RDIR = config["results_dir"] + config["run"]
CDIR = config["costs_dir"]

config.update({"git_commit": get_last_commit_message(".")})
config.update({"submodule_commit": get_last_commit_message(PYPSAEARTH_FOLDER)})

CUTOUTS_PATH = (
    "cutouts/cutout-2013-era5-tutorial.nc"
    if config["tutorial"]
    else "cutouts/cutout-2013-era5.nc"
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


if not config.get("disable_subworkflow", False):

    subworkflow pypsaearth:
        workdir:
            PYPSAEARTH_FOLDER
        snakefile:
            PYPSAEARTH_FOLDER + "/Snakefile"
        configfile:
            "./config.pypsa-earth.yaml"


if config.get("disable_subworkflow", False):

    def pypsaearth(path):
        return PYPSAEARTH_FOLDER + "/" + path


if config["enable"].get("retrieve_cost_data", True):

    rule retrieve_cost_data:
        input:
            HTTP.remote(
                f"raw.githubusercontent.com/PyPSA/technology-data/{config['costs']['version']}/outputs/costs"
                + "_{planning_horizons}.csv",
                keep_local=True,
            ),
        output:
            costs=CDIR + "costs_{planning_horizons}.csv",
        resources:
            mem_mb=5000,
        run:
            move(input[0], output[0])


rule prepare_sector_networks:
    input:
        expand(
            RDIR
            + "/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}.nc",
            **config["scenario"],
            **config["costs"]
        ),


rule override_res_all_nets:
    input:
        expand(
            RDIR
            + "/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_presec.nc",
            **config["scenario"],
            **config["costs"],
            **config["export"]
        ),


rule solve_all_networks:
    input:
        expand(
            RDIR
            + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
            **config["scenario"],
            **config["costs"],
            **config["export"]
        ),


rule prepare_ports:
    output:
        ports="data/ports.csv",  # TODO move from data to resources
    script:
        "scripts/prepare_ports.py"


rule prepare_airports:
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
        input:
            regions_onshore=pypsaearth(
                "resources/bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson"
            ),
        output:
            clustered_gas_network="resources/gas_networks/gas_network_elec_s{simpl}_{clusters}.csv",
            gas_network_fig_1="resources/gas_networks/existing_gas_pipelines_{simpl}_{clusters}.png",
            gas_network_fig_2="resources/gas_networks/clustered_gas_pipelines_{simpl}_{clusters}.png",
        script:
            "scripts/prepare_gas_network.py"


rule prepare_sector_network:
    input:
        network=RDIR
        + "/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_presec.nc",
        costs=CDIR + "costs_{planning_horizons}.csv",
        h2_cavern="data/hydrogen_salt_cavern_potentials.csv",
        nodal_energy_totals="resources/demand/heat/nodal_energy_heat_totals_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        transport="resources/demand/transport_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        avail_profile="resources/pattern_profiles/avail_profile_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        dsm_profile="resources/pattern_profiles/dsm_profile_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        nodal_transport_data="resources/demand/nodal_transport_data_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        overrides="data/override_component_attrs",
        clustered_pop_layout="resources/population_shares/pop_layout_elec_s{simpl}_{clusters}.csv",
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
        shapes_path=pypsaearth(
            "resources/bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson"
        ),
        pipelines="data_custom/pipelines.csv"
        if config["custom_data"]["gas_network"]
        else "resources/gas_networks/gas_network_elec_s{simpl}_{clusters}.csv",
    output:
        RDIR
        + "/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}.nc",
    threads: 1
    resources:
        mem_mb=2000,
    benchmark:
        (
            RDIR
            + "/benchmarks/prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}"
        )
    script:
        "scripts/prepare_sector_network.py"


rule build_ship_profile:
    output:
        ship_profile="resources/ship_profile_{h2export}TWh.csv",
    script:
        "scripts/build_ship_profile.py"


rule add_export:
    input:
        overrides="data/override_component_attrs",
        export_ports="data/export_ports.csv",
        costs=CDIR + "costs_{planning_horizons}.csv",
        ship_profile="resources/ship_profile_{h2export}TWh.csv",
        network=RDIR
        + "/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}.nc",
        shapes_path=pypsaearth(
            "resources/bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson"
        ),
    output:
        RDIR
        + "/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
    script:
        "scripts/add_export.py"


rule override_respot:
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
        network=pypsaearth("networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"),
        energy_totals="data/energy_totals_{demand}_{planning_horizons}.csv",
    output:
        RDIR
        + "/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_presec.nc",
    script:
        "scripts/override_respot.py"


rule prepare_transport_data:
    input:
        network=pypsaearth("networks/elec_s{simpl}_{clusters}.nc"),
        energy_totals_name="data/energy_totals_{demand}_{planning_horizons}.csv",
        traffic_data_KFZ="data/emobility/KFZ__count",
        traffic_data_Pkw="data/emobility/Pkw__count",
        transport_name="resources/transport_data.csv",
        clustered_pop_layout="resources/population_shares/pop_layout_elec_s{simpl}_{clusters}.csv",
        temp_air_total="resources/temperatures/temp_air_total_elec_s{simpl}_{clusters}.nc",
    output:
        # nodal_energy_totals="resources/nodal_energy_totals_s{simpl}_{clusters}.csv",
        transport="resources/demand/transport_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        avail_profile="resources/pattern_profiles/avail_profile_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        dsm_profile="resources/pattern_profiles/dsm_profile_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
        nodal_transport_data="resources/demand/nodal_transport_data_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv",
    script:
        "scripts/prepare_transport_data.py"


rule build_cop_profiles:
    input:
        temp_soil_total="resources/temperatures/temp_soil_total_elec_s{simpl}_{clusters}.nc",
        temp_soil_rural="resources/temperatures/temp_soil_rural_elec_s{simpl}_{clusters}.nc",
        temp_soil_urban="resources/temperatures/temp_soil_urban_elec_s{simpl}_{clusters}.nc",
        temp_air_total="resources/temperatures/temp_air_total_elec_s{simpl}_{clusters}.nc",
        temp_air_rural="resources/temperatures/temp_air_rural_elec_s{simpl}_{clusters}.nc",
        temp_air_urban="resources/temperatures/temp_air_urban_elec_s{simpl}_{clusters}.nc",
    output:
        cop_soil_total="resources/cops/cop_soil_total_elec_s{simpl}_{clusters}.nc",
        cop_soil_rural="resources/cops/cop_soil_rural_elec_s{simpl}_{clusters}.nc",
        cop_soil_urban="resources/cops/cop_soil_urban_elec_s{simpl}_{clusters}.nc",
        cop_air_total="resources/cops/cop_air_total_elec_s{simpl}_{clusters}.nc",
        cop_air_rural="resources/cops/cop_air_rural_elec_s{simpl}_{clusters}.nc",
        cop_air_urban="resources/cops/cop_air_urban_elec_s{simpl}_{clusters}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        "benchmarks/build_cop_profiles/s{simpl}_{clusters}"
    script:
        "scripts/build_cop_profiles.py"


rule prepare_heat_data:
    input:
        network=pypsaearth("networks/elec_s{simpl}_{clusters}.nc"),
        energy_totals_name="data/energy_totals_{demand}_{planning_horizons}.csv",
        clustered_pop_layout="resources/population_shares/pop_layout_elec_s{simpl}_{clusters}.csv",
        temp_air_total="resources/temperatures/temp_air_total_elec_s{simpl}_{clusters}.nc",
        cop_soil_total="resources/cops/cop_soil_total_elec_s{simpl}_{clusters}.nc",
        cop_air_total="resources/cops/cop_air_total_elec_s{simpl}_{clusters}.nc",
        solar_thermal_total="resources/demand/heat/solar_thermal_total_elec_s{simpl}_{clusters}.nc",
        heat_demand_total="resources/demand/heat/heat_demand_total_elec_s{simpl}_{clusters}.nc",
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
    input:
        unsd_paths="data/demand/unsd/paths/Energy_Statistics_Database.xlsx",
    output:
        energy_totals_base="data/energy_totals_base.csv",
    script:
        "scripts/build_base_energy_totals.py"


rule prepare_energy_totals:
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
    input:
        pop_layout_total="resources/population_shares/pop_layout_total.nc",
        pop_layout_urban="resources/population_shares/pop_layout_urban.nc",
        pop_layout_rural="resources/population_shares/pop_layout_rural.nc",
        regions_onshore=pypsaearth(
            "resources/bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson"
        ),
        cutout=pypsaearth(CUTOUTS_PATH),
    output:
        solar_thermal_total="resources/demand/heat/solar_thermal_total_elec_s{simpl}_{clusters}.nc",
        solar_thermal_urban="resources/demand/heat/solar_thermal_urban_elec_s{simpl}_{clusters}.nc",
        solar_thermal_rural="resources/demand/heat/solar_thermal_rural_elec_s{simpl}_{clusters}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        "benchmarks/build_solar_thermal_profiles/s{simpl}_{clusters}"
    script:
        "scripts/build_solar_thermal_profiles.py"


rule build_population_layouts:
    input:
        nuts3_shapes=pypsaearth("resources/shapes/gadm_shapes.geojson"),
        urban_percent="data/urban_percent.csv",
        cutout=pypsaearth(CUTOUTS_PATH),
    output:
        pop_layout_total="resources/population_shares/pop_layout_total.nc",
        pop_layout_urban="resources/population_shares/pop_layout_urban.nc",
        pop_layout_rural="resources/population_shares/pop_layout_rural.nc",
        gdp_layout="resources/gdp_shares/gdp_layout.nc",
    resources:
        mem_mb=20000,
    benchmark:
        "benchmarks/build_population_layouts"
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
        pop_layout_total="resources/population_shares/pop_layout_total.nc",
        pop_layout_urban="resources/population_shares/pop_layout_urban.nc",
        pop_layout_rural="resources/population_shares/pop_layout_rural.nc",
        gdp_layout="resources/gdp_shares/gdp_layout.nc",
        regions_onshore=pypsaearth(
            "resources/bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson"
        ),
        cutout=pypsaearth(CUTOUTS_PATH),
    output:
        clustered_pop_layout="resources/population_shares/pop_layout_elec_s{simpl}_{clusters}.csv",
        clustered_gdp_layout="resources/gdp_shares/gdp_layout_elec_s{simpl}_{clusters}.csv",
    resources:
        mem_mb=10000,
    benchmark:
        "benchmarks/build_clustered_population_layouts/s{simpl}_{clusters}"
    script:
        "scripts/build_clustered_population_layouts.py"


rule build_heat_demand:
    input:
        pop_layout_total="resources/population_shares/pop_layout_total.nc",
        pop_layout_urban="resources/population_shares/pop_layout_urban.nc",
        pop_layout_rural="resources/population_shares/pop_layout_rural.nc",
        regions_onshore=pypsaearth(
            "resources/bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson"
        ),
        cutout=pypsaearth(CUTOUTS_PATH),
    output:
        heat_demand_urban="resources/demand/heat/heat_demand_urban_elec_s{simpl}_{clusters}.nc",
        heat_demand_rural="resources/demand/heat/heat_demand_rural_elec_s{simpl}_{clusters}.nc",
        heat_demand_total="resources/demand/heat/heat_demand_total_elec_s{simpl}_{clusters}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        "benchmarks/build_heat_demand/s{simpl}_{clusters}"
    script:
        "scripts/build_heat_demand.py"


rule build_temperature_profiles:
    input:
        pop_layout_total="resources/population_shares/pop_layout_total.nc",
        pop_layout_urban="resources/population_shares/pop_layout_urban.nc",
        pop_layout_rural="resources/population_shares/pop_layout_rural.nc",
        regions_onshore=pypsaearth(
            "resources/bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson"
        ),
        cutout=pypsaearth(CUTOUTS_PATH),
    output:
        temp_soil_total="resources/temperatures/temp_soil_total_elec_s{simpl}_{clusters}.nc",
        temp_soil_rural="resources/temperatures/temp_soil_rural_elec_s{simpl}_{clusters}.nc",
        temp_soil_urban="resources/temperatures/temp_soil_urban_elec_s{simpl}_{clusters}.nc",
        temp_air_total="resources/temperatures/temp_air_total_elec_s{simpl}_{clusters}.nc",
        temp_air_rural="resources/temperatures/temp_air_rural_elec_s{simpl}_{clusters}.nc",
        temp_air_urban="resources/temperatures/temp_air_urban_elec_s{simpl}_{clusters}.nc",
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


if config["foresight"] == "overnight":

    rule solve_network:
        input:
            overrides="data/override_component_attrs",
            # network=RDIR
            # + "/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}.nc",
            network=RDIR
            + "/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
            costs=CDIR + "costs_{planning_horizons}.csv",
            configs=SDIR + "/configs/config.yaml",  # included to trigger copy_config rule
        output:
            RDIR
            + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
        shadow:
            "shallow"
        log:
            solver=RDIR
            + "/logs/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export_solver.log",
            python=RDIR
            + "/logs/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export_python.log",
            memory=RDIR
            + "/logs/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export_memory.log",
        threads: 25
        resources:
            mem_mb=config["solving"]["mem"],
        benchmark:
            (
                RDIR
                + "/benchmarks/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export"
            )
        script:
            "scripts/solve_network.py"


rule make_summary:
    input:
        overrides="data/override_component_attrs",
        networks=expand(
            RDIR
            + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
            **config["scenario"],
            **config["costs"],
            **config["export"]
        ),
        costs=CDIR + "costs_{}.csv".format(config["scenario"]["planning_horizons"][0]),
        plots=expand(
            RDIR
            + "/maps/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}-costs-all_{planning_horizons}_{discountrate}_{demand}_{h2export}export.pdf",
            **config["scenario"],
            **config["costs"],
            **config["export"]
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
        + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
    output:
        map=RDIR
        + "/maps/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}-costs-all_{planning_horizons}_{discountrate}_{demand}_{h2export}export.pdf",
    threads: 2
    resources:
        mem_mb=10000,
    benchmark:
        (
            RDIR
            + "/benchmarks/plot_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export"
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


rule build_industrial_database:
    output:
        industrial_database="data/industrial_database.csv",
    script:
        "scripts/build_industrial_database.py"


rule prepare_db:
    input:
        network=RDIR
        + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
    output:
        db=RDIR
        + "/summaries/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}-costs-all_{planning_horizons}_{discountrate}_{demand}_{h2export}export.csv",
    threads: 2
    resources:
        mem_mb=10000,
    benchmark:
        (
            RDIR
            + "/benchmarks/prepare_db/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export"
        )
    script:
        "scripts/prepare_db.py"


rule run_test:
    params:
        dummy="This is a dummy parameter to satisfy snakefmt",
    run:
        import yaml

        with open(PYPSAEARTH_FOLDER + "/config.tutorial.yaml") as file:
            config_pypsaearth = yaml.full_load(file)
            config_pypsaearth["retrieve_databundle"] = {"show_progress": False}
            config_pypsaearth["electricity"]["extendable_carriers"]["Store"] = []
            config_pypsaearth["electricity"]["extendable_carriers"]["Link"] = []
            config_pypsaearth["electricity"]["co2limit"] = 7.75e7

            with open("./config.pypsa-earth.yaml", "w") as wfile:
                yaml.dump(config_pypsaearth, wfile)

        shell("cp test/config.test1.yaml config.yaml")
        shell("snakemake --cores all solve_all_networks --forceall")



rule clean:
    run:
        shell("rm -r " + PYPSAEARTH_FOLDER + "/resources")
        shell("rm -r " + PYPSAEARTH_FOLDER + "/networks")


rule build_industrial_distribution_key:  #default data
    input:
        regions_onshore=pypsaearth(
            "resources/bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson"
        ),
        clustered_pop_layout="resources/population_shares/pop_layout_elec_s{simpl}_{clusters}.csv",
        clustered_gdp_layout="resources/gdp_shares/gdp_layout_elec_s{simpl}_{clusters}.csv",
        industrial_database="data/industrial_database.csv",
        shapes_path=pypsaearth(
            "resources/bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson"
        ),
    output:
        industrial_distribution_key="resources/demand/industrial_distribution_key_elec_s{simpl}_{clusters}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/build_industrial_distribution_key_elec_s{simpl}_{clusters}"
    script:
        "scripts/build_industrial_distribution_key.py"


rule build_base_industry_totals:  #default data
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
    input:
        industrial_distribution_key="resources/demand/industrial_distribution_key_elec_s{simpl}_{clusters}.csv",
        #industrial_production_per_country_tomorrow="resources/demand/industrial_production_per_country_tomorrow_{planning_horizons}_{demand}.csv",
        #industrial_production_per_country="data/industrial_production_per_country.csv",
        base_industry_totals="resources/demand/base_industry_totals_{planning_horizons}_{demand}.csv",
        industrial_database="data/industrial_database.csv",
        costs=CDIR + "costs_{}.csv".format(config["scenario"]["planning_horizons"][0]),
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


if config["enable"].get("retrieve_irena", True):

    rule retrieve_irena:
        output:
            offwind="data/existing_infrastructure/offwind_capacity_IRENA.csv",
            onwind="data/existing_infrastructure/onwind_capacity_IRENA.csv",
            solar="data/existing_infrastructure/solar_capacity_IRENA.csv",
        # log:
        #     logs("retrieve_irena.log"),
        resources:
            mem_mb=1000,
        script:
            "./scripts/retrieve_irena.py"


if config["foresight"] == "myopic":

    rule add_existing_baseyear:
        params:
            baseyear=config["scenario"]["planning_horizons"][0],
            sector=config["sector"],
            existing_capacities=config["existing_capacities"],
            costs=config["costs"],
        input:
            network=RDIR
            + "/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
            powerplants=pypsaearth("resources/powerplants.csv"),
            busmap_s=pypsaearth("resources/bus_regions/busmap_elec_s{simpl}.csv"),
            busmap=pypsaearth(
                "resources/bus_regions/busmap_elec_s{simpl}_{clusters}.csv"
            ),
            clustered_pop_layout="resources/population_shares/pop_layout_elec_s{simpl}_{clusters}.csv",
            costs=CDIR
            + "costs_{}.csv".format(config["scenario"]["planning_horizons"][0]),
            cop_soil_total="resources/cops/cop_soil_total_elec_s{simpl}_{clusters}.nc",
            cop_air_total="resources/cops/cop_air_total_elec_s{simpl}_{clusters}.nc",
            existing_heating_distribution="data/existing_infrastructure/existing_heating_raw.csv",
            existing_solar="data/existing_infrastructure/solar_capacity_IRENA.csv",
            existing_onwind="data/existing_infrastructure/onwind_capacity_IRENA.csv",
            existing_offwind="data/existing_infrastructure/offwind_capacity_IRENA.csv",
        output:
            RDIR
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
            RDIR
            + "/logs/add_existing_baseyear_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.log",
        benchmark:
            (
                RDIR
                + "/benchmarks/add_existing_baseyear/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export"
            )
        script:
            "scripts/add_existing_baseyear.py"

    def input_profile_tech_brownfield(w):
        return {
            f"profile_{tech}": pypsaearth(
                f"resources/renewable_profiles/profile_{tech}.nc"
            )
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
            unpack(input_profile_tech_brownfield),
            simplify_busmap=pypsaearth("resources/bus_regions/busmap_elec_s{simpl}.csv"),
            cluster_busmap=pypsaearth(
                "resources/bus_regions/busmap_elec_s{simpl}_{clusters}.csv"
            ),
            network=RDIR
            + "/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
            network_p=solved_previous_horizon,  #solved network at previous time step
            costs=CDIR + "costs_{planning_horizons}.csv",
            cop_soil_total="resources/cops/cop_soil_total_elec_s{simpl}_{clusters}.nc",
            cop_air_total="resources/cops/cop_air_total_elec_s{simpl}_{clusters}.nc",
        output:
            RDIR
            + "/prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
        threads: 4
        resources:
            mem_mb=10000,
        log:
            RDIR
            + "/logs/add_brownfield_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.log",
        benchmark:
            (
                RDIR
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
            network=RDIR
            + "/prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
            costs=CDIR + "costs_{planning_horizons}.csv",
            configs=SDIR + "/configs/config.yaml",  # included to trigger copy_config rule
        output:
            network=RDIR
            + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
            # config=RDIR
            # + "/configs/config.elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.yaml",
        shadow:
            "shallow"
        log:
            solver=RDIR
            + "/logs/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export_solver.log",
            python=RDIR
            + "/logs/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export_python.log",
            memory=RDIR
            + "/logs/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export_memory.log",
        threads: 25
        resources:
            mem_mb=config["solving"]["mem"],
        benchmark:
            (
                RDIR
                + "/benchmarks/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export"
            )
        script:
            "./scripts/solve_network.py"

    rule solve_all_networks_myopic:
        input:
            networks=expand(
                RDIR
                + "/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc",
                **config["scenario"],
                **config["costs"],
                **config["export"],
            ),
