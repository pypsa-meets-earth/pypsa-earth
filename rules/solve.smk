# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later


if config["monte_carlo"]["options"].get("add_to_snakefile", False) != True:

    rule solve_network:
        params:
            solving=config["solving"],
            augmented_line_connection=config["augmented_line_connection"],
            policy_config=config["policy_config"],
        input:
            network=rules.prepare_network.output.network,
            agg_p_nom_minmax=config["electricity"]["agg_p_nom_limits"]["file"],  # ensure the CSV with capacity constraints is copied into the shadow directory (needed on Windows, since shadowed scripts can’t access files outside `input`)
        output:
            network="results/"
            + RDIR
            + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        log:
            solver=os.path.normpath(
                "logs/"
                + RDIR
                + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_solver.log"
            ),
            python="logs/"
            + RDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_python.log",
        benchmark:
            (
                "benchmarks/"
                + RDIR
                + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
            )
        threads: 20
        resources:
            mem=memory,
        shadow:
            "copy-minimal" if os.name == "nt" else "shallow"
        script:
            scripts("solve_network.py")


if config["monte_carlo"]["options"].get("add_to_snakefile", False) == True:

    rule solve_monte:
        input:
            expand(rules.monte_carlo.output.network, **config["scenario"]),

    rule solve_network:
        params:
            solving=config["solving"],
            augmented_line_connection=config["augmented_line_connection"],
            policy_config=config["policy_config"],
        input:
            network=rules.monte_carlo.output.network,
            agg_p_nom_minmax=config["electricity"]["agg_p_nom_limits"]["file"],  # ensure the CSV with capacity constraints is copied into the shadow directory (needed on Windows, since shadowed scripts can’t access files outside `input`)
        output:
            network="results/"
            + RDIR
            + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}.nc",
        log:
            solver=os.path.normpath(
                "logs/"
                + RDIR
                + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}_solver.log"
            ),
            python="logs/"
            + RDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}_python.log",
            memory="logs/"
            + RDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}_memory.log",
        benchmark:
            (
                "benchmarks/"
                + RDIR
                + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}"
            )
        threads: 20
        resources:
            mem_mb=memory,
        shadow:
            "copy-minimal" if os.name == "nt" else "shallow"
        script:
            scripts("solve_network.py")

    rule solve_all_networks_monte:
        input:
            expand(rules.solve_network.output.network, **config["scenario"]),


if config["foresight"] == "overnight":

    rule solve_sector_network:
        params:
            solving=config["solving"],
            augmented_line_connection=config["augmented_line_connection"],
            policy_config=config["policy_config"],
        input:
            network=rules.add_export.output.network,
            costs=rules.process_cost_data.output.costs.format(
                year="{planning_horizons}", scope="sec"
            ),
            configs=rules.copy_config.output.config,  # included to trigger copy_config rule
            agg_p_nom_minmax=config["electricity"]["agg_p_nom_limits"]["file"],  # ensure the CSV with capacity constraints is copied into the shadow directory (needed on Windows, since shadowed scripts can’t access files outside `input`)
        output:
            network=RESDIR
            + "postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}.nc",
        shadow:
            "copy-minimal" if os.name == "nt" else "shallow"
        log:
            solver="logs/"
            + SECDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_solver.log",
            python="logs/"
            + SECDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_python.log",
            memory="logs/"
            + SECDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_memory.log",
        threads: 25
        resources:
            mem_mb=config["solving"]["mem"],
        benchmark:
            (
                RESDIR
                + "benchmarks/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}"
            )
        script:
            scripts("solve_network.py")


if config["foresight"] == "myopic":

    HEAT_BASEYEAR = {
        "cop_soil_total": rules.build_cop_profiles.output.cop_soil_total,
        "cop_air_total": rules.build_cop_profiles.output.cop_air_total,
        "existing_heating_distribution": rules.build_existing_heating_distribution.output.existing_heating_distribution,
    }

    rule add_existing_baseyear:
        params:
            baseyear=config["scenario"]["planning_horizons"][0],
            sector=config["sector"],
            existing_capacities=config["existing_capacities"],
            costs=config["costs"],
        input:
            **branch(sector_enable["heat"], HEAT_BASEYEAR),
            network=rules.add_export.output.network,
            powerplants=rules.build_powerplants.output.powerplants,
            busmap_s=rules.simplify_network.output.busmap,
            busmap=rules.cluster_network.output.busmap,
            # clustered_pop_layout="resources/"
            # + SECDIR
            # + "population_shares/pop_layout_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
            costs=rules.process_cost_data.output.costs.format(
                year="{planning_horizons}", scope="sec"
            ),
        output:
            network=RESDIR
            + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}.nc",
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
            + "logs/add_existing_baseyear_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}.log",
        benchmark:
            RESDIR
            +"benchmarks/add_existing_baseyear/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}"
        script:
            scripts("add_existing_baseyear.py")

    def input_profile_tech_brownfield(w):
        return {
            f"profile_{tech}": rules.build_renewable_profiles.output.profile.format(
                technology=tech
            )
            for tech in config["electricity"]["renewable_carriers"]
            if tech != "hydro"
        }

    def solved_previous_horizon(w):
        planning_horizons = config["scenario"]["planning_horizons"]
        i = planning_horizons.index(int(w.planning_horizons))
        prev_planning_horizon = str(planning_horizons[i - 1])
        return rules.solve_network_myopic.output.network.format(
            **{**dict(w), "planning_horizons": prev_planning_horizon}
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
            simplify_busmap=rules.simplify_network.output.busmap,
            cluster_busmap=rules.cluster_network.output.busmap,
            network=rules.add_export.output.network,
            network_p=solved_previous_horizon,  #solved network at previous time step
            costs=rules.process_cost_data.output.costs.format(
                year="{planning_horizons}", scope="sec"
            ),
            cop_soil_total=rules.build_cop_profiles.output.cop_soil_total,
            cop_air_total=rules.build_cop_profiles.output.cop_air_total,
        output:
            network=RESDIR
            + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}.nc",
        threads: 4
        resources:
            mem_mb=10000,
        log:
            RESDIR
            + "logs/add_brownfield_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}.log",
        benchmark:
            (
                RESDIR
                + "benchmarks/add_brownfield/elec_s{simpl}_ec_{clusters}_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}"
            )
        script:
            scripts("add_brownfield.py")

    ruleorder: add_existing_baseyear > add_brownfield

    rule solve_network_myopic:
        params:
            solving=config["solving"],
            foresight=config["foresight"],
            planning_horizons=config["scenario"]["planning_horizons"],
            co2_sequestration_potential=config["scenario"].get(
                "co2_sequestration_potential", 200
            ),
            augmented_line_connection=config["augmented_line_connection"],
            policy_config=config["policy_config"],
        input:
            network=RESDIR
            + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}.nc",
            costs=rules.process_cost_data.output.costs.format(
                year="{planning_horizons}", scope="sec"
            ),
            configs=rules.copy_config.output.config,  # included to trigger copy_config rule
            agg_p_nom_minmax=config["electricity"]["agg_p_nom_limits"]["file"],  # ensure the CSV with capacity constraints is copied into the shadow directory (needed on Windows, since shadowed scripts can’t access files outside `input`)
        output:
            network=RESDIR
            + "postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}.nc",
            # config=RESDIR
            # + "configs/config.elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}.yaml",
        shadow:
            "copy-minimal" if os.name == "nt" else "shallow"
        log:
            solver="logs/"
            + SECDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_solver.log",
            python="logs/"
            + SECDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_python.log",
            memory="logs/"
            + SECDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_memory.log",
        threads: 25
        resources:
            mem_mb=config["solving"]["mem"],
        benchmark:
            (
                RESDIR
                + "benchmarks/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}"
            )
        script:
            scripts("solve_network.py")


rule run_scenario:
    input:
        diff_config="configs/scenarios/config.{scenario_name}.yaml",
    output:
        touchfile=touch("results/{scenario_name}/scenario.done"),
        copyconfig="results/{scenario_name}/config.yaml",
    threads: 1
    resources:
        mem_mb=5000,
    run:
        from subprocess import run

        import yaml
        from build_test_configs import create_test_config

        # get base configuration file from diff config
        with open(input.diff_config) as f:
            base_config_path = (
                yaml.full_load(f)
                .get("run", {})
                .get("base_config", "config.default.yaml")
            )

            # Ensure the scenario name matches the name of the configuration
        create_test_config(
            input.diff_config,
            {"run": {"name": wildcards.scenario_name}},
            input.diff_config,
        )
        # merge the default config file with the difference
        create_test_config(base_config_path, input.diff_config, "config.yaml")
        run(
            "snakemake -j all solve_all_networks --rerun-incomplete",
            shell=True,
            check=not config["run"]["allow_scenario_failure"],
        )
        run(
            "snakemake -j1 make_statistics --force",
            shell=True,
            check=not config["run"]["allow_scenario_failure"],
        )
        copyfile("config.yaml", output.copyconfig)



rule run_all_scenarios:
    input:
        expand(
            rules.run_scenario.output.touchfile,
            scenario_name=[
                c.stem.replace("config.", "")
                for c in Path("configs/scenarios").glob("config.*.yaml")
            ],
        ),
