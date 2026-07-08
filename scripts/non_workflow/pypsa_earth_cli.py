# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
This script provides the basis to run a command line interface (CLI) to help users navigate through PyPSA-Earth

The CLI has the following modules:
1. Tutorial - A tutorial based module developed in sync with use-case documentation
2. Edit config - Feature to edit config parameters without direct exposure to the same
3. Retrieve databundles - Feature to independently retrieve databundles required for a PyPSA-Earth model run
4. Run snakemake - Feature to run a snakemake workflow

"""
import os
import subprocess
import sys
from enum import Enum

import typer
import yaml
from InquirerPy import get_style, inquirer
from InquirerPy.base import Choice
from rich.columns import Columns
from rich.console import Console
from rich.panel import Panel

app = typer.Typer(help="CLI to change config entries in PyPSA-Earth and run the model")
console = Console()

# Use get_style to create a proper InquirerPyStyle object
custom_menu_style = get_style(
    {
        "questionmark": "hidden",
        "question": "yellow",  # Title color inside the top border frame
        "pointer": "#00ffff bold",
        "choice": "#808080",
        "selected": "#ff00ff bold underline",
        "frame": "bold white",  # Border lines
    },
    style_override=False,
)


def display_choice_menu(title: str, options: list, required_height: int) -> str:
    """
    Display a menu

    Parameters
    ----------
    title: str
        Title for the menu
    options: list
        List of options to be displayed
    required_height: int
        Height of the border box for the option list

    Returns
    -------
    choice: str
        Choice selected by the user
    """

    choice = inquirer.select(
        message=title,
        choices=options,
        pointer="👉",
        style=custom_menu_style,
        border=True,
        max_height=required_height,
    ).execute()

    return choice


def ask(text: str, default: str = "") -> str:
    """
    Styling for the typer prompt

    Parameters
    ----------
    text: str
        Text to be shown in the prompt
    default: str
        The default values for the prompt

    Returns:
    str
        Stylized text for the prompt
    """
    styled_text = typer.style(f"➔ {text}", fg=typer.colors.YELLOW, bold=True)
    if default != "":
        return typer.prompt(styled_text, default=default)
    else:
        return typer.prompt(styled_text)


def exit_message() -> None:
    console.print(style="dim")
    console.print("[bold magenta] 👋 Exiting the application. [/bold magenta]")
    raise typer.Exit()


def load_config_file(config_path: str) -> dict:
    """
    Load configuration file

    Parameters
    ----------
    config_path: str
        Path to the config file to be loaded

    Returns
    -------
    dict:
        loaded config data
    """
    with open(config_path, "r") as file:
        try:
            config = yaml.safe_load(file)
            return config
        except yaml.YAMLError as error:
            print(f"Error parsing YAML file: {error}")


def save_config_file(config_path: str, config_data: dict[dict]) -> None:
    """
    Save a nested dictionary as a config file

    Parameters
    ----------
    config_path: str
        Path to save the config file to
    config_data: dict[dict]
        Nested dictionary of updated config values

    Returns
    -------
    None
    """
    with open(config_path, "w") as file:
        yaml.safe_dump(config_data, file, default_flow_style=False)
        console.print(f"[bold green] ✔ Saved the file to {config_path} [/bold green]")
        console.print(style="dim")


def flatten_dict(
    nested_dict: dict, separator: str = ".", current_path: str = ""
) -> dict:
    """
    Recursively flattens a nested dictionary using a separator.

    Parameters
    ----------
    nested_dict: dict
        Nested dictionary of config values
    separator: str
        Separator when flattening a dict

    Returns
    -------
    dict
        Flattening of nested dictionary
    """
    flat_map = {}

    for key, value in nested_dict.items():
        # Build the new key path string
        new_path = f"{current_path}{separator}{key}" if current_path else str(key)

        # If the value is another dictionary, dig deeper recursively
        if isinstance(value, dict) and value:
            flat_map.update(flatten_dict(value, separator, new_path))
        else:
            # Leaf node reached, map it to the flattened key
            flat_map[new_path] = value

    return flat_map


def unflatten_dict(flat_dict: dict, separator: str = ".") -> dict:
    """
    Converts a flat dictionary with delimited keys back into a nested dictionary.

    Parameters
    ----------
    flat_dict: dict
        Flattened dictionary
    separator: str
        Separator string

    Returns
    --------
    Unflattened nested dictionary
    """
    nested_dict = {}

    for flat_key, value in flat_dict.items():
        # Split the path string into individual keys
        keys = flat_key.split(separator)

        # Traverse down or build the nested dictionary layers
        current_layer = nested_dict
        for key in keys[:-1]:
            # If the layer doesn't exist yet, initialize an empty dict
            if key not in current_layer or not isinstance(current_layer[key], dict):
                current_layer[key] = {}
            current_layer = current_layer[key]
        # Assign the actual value to the final leaf key
        current_layer[keys[-1]] = value

    return nested_dict


def display_user_groups() -> tuple:
    """
    Menu to display user groups to the user from those defined in `user_groups_cli.yaml`.

    Returns
    -------
    user_groups_config: dict[str, list[str]]
        Dictionary of config parameters associated with each user group
    user_group: str
        Choice of user group
    """
    # Load user groups
    user_groups_path = "configs/pypsa_earth_cli/user_groups.yaml"
    user_groups_config = load_config_file(user_groups_path)

    # Prompt user to identify his/her user group
    options = []
    for index, (key, value) in enumerate(user_groups_config.items(), start=1):
        display_label = f"{key} : {value}"
        options.append(Choice(value=key, name=display_label))

    # To add more config parameters to a user group, modify user_groups.yaml
    title = "Select user group"
    user_group = display_choice_menu(title, options, len(options) + 1)
    return user_groups_config, user_group


# Currently not in use
def display_config_files() -> dict[dict]:
    """
    Select and load config file data

    Returns
    -------
    dict[dict]
        Nested config parameters
    """
    config_files_list = [
        f for f in os.listdir(".") if f.endswith(".yaml") and f.startswith("config")
    ]

    config_path = display_choice_menu(
        "Select config file", config_files_list, len(config_files_list) + 1
    )

    config = load_config_file(config_path)
    return config


@app.command("config-setup")
def config_setup():
    user_groups_config, user_group = display_user_groups()

    # Load config file
    config_path = "config.default.yaml"
    config = load_config_file(config_path)

    # Filter config file for the required params based on user group
    config_options = user_groups_config[user_group]

    # Fetch required config file parameters based on user group
    reqd_config = {x: config[x] for x in config_options}

    # Create a copy of the required config to store updated values
    updated_config = reqd_config.copy()

    config_options.append("return to main menu")
    return_to_main = False
    choice = ""
    while not return_to_main:
        # Use case choice
        required_height = len(config_options) + 2

        # Iterate through parent config
        title = "🗺️  Select config option to modify"
        choice = display_choice_menu(title, config_options, required_height)

        if choice == "return to main menu":
            return_to_main = True
            console.print("[bold blue]⏳ Returning to main menu [/bold blue]")
            break

        # Iterate through the child config parameters, e.g., for the `scenario` config option, the child params are `simpl`,`ll`,`clusters` etc.
        if isinstance(reqd_config[choice], dict):
            return_to_config_setup = False
            while not return_to_config_setup:
                flattened_options = flatten_dict(reqd_config[choice])

                subchoice_options = []
                for index, (key, value) in enumerate(
                    flattened_options.items(), start=1
                ):
                    # The 'name' controls what the user sees in the box.
                    # The 'value' is what gets returned to Python upon selection.
                    display_label = f"{index}. {key} : {value}"
                    subchoice_options.append(Choice(value=key, name=display_label))

                subchoice_options.append(
                    Choice(
                        value="return",
                        name=f"{index+1}. Go back to parent config options",
                    )
                )
                required_height = len(subchoice_options) + 2
                title = "🗺️  Select config option to modify"
                subchoice = display_choice_menu(
                    title, subchoice_options, required_height
                )

                if subchoice == "return":
                    return_to_config_setup = True
                    continue

                updated_value = ask(
                    f"Enter the value of {subchoice} to update in the config file"
                )
                if isinstance(updated_config[choice][subchoice], list):
                    updated_value = updated_value.split(",")
                updated_config[choice][subchoice] = updated_value
        else:
            # If no nested params exist for a particular config option, directly ask for the value to be updated. E.g., countries
            updated_value = ask(
                f"Enter the value of {choice} to update in the config file"
            )
            if isinstance(updated_config[choice], list):
                # Splitting by separator to allow for multiple values to be entered for a config option. E.g., countries
                updated_value = updated_value.split(",")
            updated_config[choice] = updated_value

    # Save updated config file
    config_save_path = "config.cli_updated.yaml"
    save_config_file(config_save_path, unflatten_dict(updated_config))
    display_main_menu()


def display_main_menu() -> None:
    """
    Interactive menu to be displayed as a starting point for the CLI
    """
    menu_items = [
        {
            "num": "1",
            "name": "Tutorial",
            "desc": "Use-case guide for PyPSA-Kazakhstan",
        },
        {
            "num": "2",
            "name": "Configuration Setup",
            "desc": "Customize configuration parameters for PyPSA-Earth model run",
        },
        {
            "num": "3",
            "name": "Data Retrieval",
            "desc": "Fetch external data required for PyPSA-Earth model run",
        },
        {"num": "4", "name": "Run Model", "desc": "Run PyPSA-Earth model"},
        {"num": "5", "name": "Exit", "desc": "Close the application"},
    ]

    panels = []
    for item in menu_items:
        panel_content = f"[bold cyan]{item['num']}[/bold cyan] | [bold white]{item['name']}[/bold white]\n[dim]{item['desc']}[/dim]"
        panels.append(Panel(panel_content, expand=False, border_style="green"))

    # Print the menu title and options
    console.rule("[bold magenta]📊 MAIN MENU [/bold magenta]")
    console.print(Columns(panels, padding=(1, 2)))
    choice = ask("Select option 1-5 to proceed further")

    if choice == "2":
        config_setup()
    elif choice == "3":
        subprocess.run([sys.executable, "scripts/non_workflow/databundle_cli.py"])
        display_main_menu()
    elif choice == "5":
        exit_message()
    else:
        console.print(
            "[bold magenta] Feature still under development. Please check again later [/bold magenta]"
        )
        exit_message()


@app.callback(invoke_without_command=True)
def main(ctx: typer.Context):

    if ctx.invoked_subcommand is not None:
        return

    # Display the welcome message
    console.rule(
        "[bold magenta]🚀 Welcome to PyPSA-Earth config CLI application! [/bold magenta]"
    )
    console.print(
        "[white]PyPSA-Earth is the first open-source global cross-sectoral energy system model "
        "with high spatial and temporal resolution. The workflow provide capabilities for modelling the energy systems of"
        " any country in the world, enabling large-scale collaboration and transparent analysis for an inclusive and "
        "sustainable energy future. PyPSA-Earth is suitable for both operational studies and capacity expansion studies. "
        "Its sector-coupled modeling capabilities enable features for the detailed optimization of multi-energy systems, "
        "covering electricity, heating, transport, industry, hydrogen and more. [/white]"
    )
    console.print(style="dim")

    display_main_menu()


if __name__ == "__main__":
    app()
