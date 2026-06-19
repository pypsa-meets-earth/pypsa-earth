from enum import Enum
from pathlib import Path

import typer
import yaml
from rich.columns import Columns
from rich.console import Console
from rich.panel import Panel

app = typer.Typer(help="CLI to change config entries in PyPSA-Earth and run the model")
console = Console()


class UserGroup(str, Enum):
    beg = "Beginner"
    mod = "Modeller"
    client = "Client"
    custom = "Custom"


def ask(text: str, default: str = "") -> str:
    styled_text = typer.style(f"➔ {text}", fg=typer.colors.YELLOW, bold=True)
    if default != "":
        return typer.prompt(styled_text, default=default)
    else:
        return typer.prompt(styled_text)


def load_config_file(config_path: str) -> dict:
    with open(config_path, "r") as file:
        try:
            config = yaml.safe_load(file)
            return config
        except yaml.YAMLError as error:
            print(f"Error parsing YAML file: {error}")


def save_config_file(config_path: str, config_data: dict[dict]) -> None:
    with open(config_path, "w") as file:
        yaml.safe_dump(config_data, file, default_flow_style=False)
        print(f"Saved the file to {config_path}")


def flatten_dict(
    nested_dict: dict, separator: str = ".", current_path: str = ""
) -> dict:
    """Recursively flattens a nested dictionary using a path separator."""
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
    """Converts a flat dictionary with delimited keys back into a nested dictionary."""
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


@app.command("config-setup")
def config_setup():
    console.print(style="dim")
    console.rule("[bold magenta] EDIT CONFIG PARAMETERS [/bold magenta]")

    # Load user groups
    user_groups_path = "user_groups.yaml"
    user_groups_config = load_config_file(user_groups_path)
    # Prompt user to identify his/her user group
    user_group = ask("Select user group", default=[x.value for x in UserGroup])

    # Prompt user for config file locations
    config_path = ask("Enter config file name to modify", default="config.default.yaml")
    console.print(style="dim")
    try:
        config = load_config_file(config_path)
    except FileNotFoundError:
        console.print(
            "[bold red] No such file found. \n Exitting the application [/bold red]"
        )
        raise typer.Exit()

    # Filter config file for the required params based on user group
    config_options = user_groups_config[user_group]
    # Read config file
    try:
        reqd_config = {x: config[x] for x in config_options}
    except KeyError:
        console.print(
            "[bold red] Mismatch in choice of user group and config file. \n Exitting the application [/bold red]"
        )
        raise typer.Exit()

    console.print("[bold]Instructions[/bold]\n")
    console.print(
        "If a config entry is a nested dictionary, the child param name is represented the parent param name separated by a '.'.\n"
    )
    console.print(style="dim")
    console.print("For example, \nenable:\n\tretrieve_databundle:\n")
    console.print(style="dim")
    console.print("is represented as [bold]enable.retrieve_databundle[/bold].")
    console.print(style="dim")
    console.print(style="dim")

    # Flatten nested dictionaries in the config parameters
    flattened_options = flatten_dict(reqd_config)

    # Iterate through the config parameters
    updated_config = flattened_options.copy()
    for option in flattened_options:
        updated_value = ask(
            f"Enter new value for {option} [Press Enter to skip]",
            default=flattened_options[option],
        )
        console.print(style="dim")

        if updated_value != flattened_options[option]:
            updated_config[option] = updated_value

    # Save updated config file
    config_save_path = ask(
        f"Enter name for updated config file [Press Enter to skip]",
        default="config.cli_updated.yaml",
    )
    save_config_file(config_save_path, unflatten_dict(updated_config))


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

    menu_items = [
        {
            "num": "1",
            "name": "Setup & activate environment",
            "desc": "Setup python environment using conda",
        },
        {
            "num": "2",
            "name": "Edit config parameters",
            "desc": "Edit config parameters",
        },
        {
            "num": "3",
            "name": "Run model",
            "desc": "Trigger snakemake workflow to run PyPSA-Earth model",
        },
        {"num": "4", "name": "Exit", "desc": "Close the application"},
    ]

    panels = []
    for item in menu_items:
        panel_content = f"[bold cyan]{item['num']}[/bold cyan] | [bold white]{item['name']}[/bold white]\n[dim]{item['desc']}[/dim]"
        panels.append(Panel(panel_content, expand=False, border_style="green"))

    # Print the menu title and options
    console.rule("[bold magenta] MAIN MENU [/bold magenta]")
    console.print(Columns(panels, padding=(1, 2)))
    choice = ask("Select option 1-4 to proceed further")

    if choice == "2":
        config_setup()
    else:
        console.print(
            "[bold magenta] Feature still under development. Please check again later [/bold magenta]"
        )
        raise typer.Exit()


if __name__ == "__main__":
    app()
