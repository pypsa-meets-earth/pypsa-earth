# SPDX-FileCopyrightText: Open Energy Transition gGmbH and contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
import os
import glob
import shutil

import yaml

# ---------------------- Utility Functions ----------------------


def deep_update(original, updates):
    if not updates:
        return

    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(original.get(key), dict):
            deep_update(original[key], value)
        else:
            original[key] = value


def select_scenario(scenarios, name="scenarios", info=True):
    print(f"Available {name}:")
    for i, scenario in enumerate(scenarios, 1):
        print(f"{i}. {scenario}")

    if info:
        print(f"(It must be placed within the config/ folder, has the letter '{name}' in the beginning and '.yaml' in the end)")
    print("(You can also select the list number)")

    while True:
        answer = input(f"Select a {name} (or press Enter to cancel): ").strip()
        if not answer:
            return None
        if answer in scenarios:
            print(f"'{answer}' is selected.")
            return answer
        elif answer.isdigit():
            scenario_selected = scenarios[int(answer) - 1]
            print(f"'{scenario_selected}' is selected.")
            return scenario_selected
        else:
            print(f"'{answer}' not found in the list. Please try again.")


def select_profile():
    while True:
        answer = input("Is this in a computer cluster [y/n]?: ").strip().lower()
        if not answer:
            return None
        elif answer in ["y", "yes"]:
            return "--profile slurm"
        elif answer in ["n", "no"]:
            cpu = (
                input("How many CPUs do you want to use (all or a number)? ")
                .strip()
                .lower()
            )
            if cpu == "all":
                return "-call"
            elif cpu.isdigit():
                return f"-c{cpu}"
        print(f"'{answer}' is not a valid option. Please try again.")


def select_multiruns():
    while True:
        answer = (
            input("is there more runs that you want to make [y/n]?: ").strip().lower()
        )
        if not answer:
            return None
        elif answer in ["y", "yes"]:
            return True
        elif answer in ["n", "no"]:
            return False
        print(f"'{answer}' is not a valid option. Please try again.")


# ---------------------- Main Script ----------------------


def main():
    base_config_yaml = glob.glob("configs/config*.yaml")

    selected_base_config_yaml = select_scenario(base_config_yaml, name="config")
    if not selected_base_config_yaml:
        print("\nOperation cancelled by user.")
        return

    print("\n=================================================================")

    scenario_yaml = glob.glob("configs/scenarios*.yaml")

    selected_scenario_yaml = select_scenario(scenario_yaml, name="scenarios")
    if not selected_scenario_yaml:
        print("\nOperation cancelled by user.")
        return

    print("\n=================================================================")
    # Load scenario configuration
    with open(selected_scenario_yaml) as file:
        config_s = yaml.safe_load(file)

    # Extract scenarios with '--'
    scenario_list = list(config_s.keys())

    scenario_pair = []
    while True:
        # User input
        selected_scenario = select_scenario(scenario_list, name="scenario", info=False)
        if not selected_scenario:
            print("\nOperation cancelled by user.")
            return

        scenario_pair += [selected_scenario]

        print("\nFinal selection:")
        print(f"base_config: {selected_base_config_yaml}")
        print(f"scenarios_config: {selected_scenario_yaml}")
        for selected_scenario in scenario_pair:
            print(f"- {selected_scenario}")

        selected_multiruns = select_multiruns()
        if not selected_multiruns:
            break

    print("\n=================================================================")
    selected_profile = select_profile()
    if not selected_profile:
        print("\nOperation cancelled by user.")
        return

    print(f"Profile: {selected_profile}")

    # Iterate runs
    for selected_scenario in scenario_pair:
        print("\n=================================================================")
        print("Currently running:")
        print(f"Baseline: {selected_scenario}")

        with open(selected_base_config_yaml) as file:
            config = yaml.safe_load(file)

        deep_update(config, config_s[selected_scenario])

        config_update = {
            "run": {
                "name": selected_scenario, 
                "sector_name": selected_scenario
            }
        }

        deep_update(config, config_update)

        # Write temp config
        with open("config.temp.yaml", "w") as file:
            yaml.safe_dump(config, file, default_flow_style=False)

        run_cmd = f"snakemake {selected_profile} solve_sector_networks --configfile config.temp.yaml"
        os.system(run_cmd)

        temp_file = "config.temp.yaml"

        if os.path.exists(temp_file):
            os.remove(temp_file)
            print(f"Deleted file: {temp_file}")
        else:
            print(f"File does not exist: {temp_file}")


if __name__ == "__main__":
    main()
