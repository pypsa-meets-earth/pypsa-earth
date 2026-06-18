import typer
from rich.console import Console
from rich.panel import Panel
from pathlib import Path
from typing_extensions import Annotated
import yaml
from enum import Enum

app = typer.Typer(help="CLI to change config entries in PyPSA-Earth and run the model")
console=Console()

class UserGroup(str, Enum):
    beg = "Beginner"
    mod = "Modeller"
    client = "Client"
    custom = "Custom"

def load_config_file(config_path):
    with open(config_path, "r") as file:
        try:
            config = yaml.safe_load(file)
            return config
        except yaml.YAMLError as error:
            print(f"Error parsing YAML file: {error}")


def save_config_file(config_path, config_data):
    with open(config_path, "w") as file:
        yaml.safe_dump(config_data, file, default_flow_style=False)
        print(f"Saved the file to {config_path}")



def flatten_dict(nested_dict: dict, separator: str = ".", current_path: str = "") -> dict:
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

@app.command()
def config_setup(user_group:Annotated[UserGroup,typer.Option(prompt="Select user group")]):
    config_path="config.default.yaml"
    config=load_config_file(config_path)

    user_groups_path="user_groups.yaml"
    user_groups_config=load_config_file(user_groups_path)

    config_options=user_groups_config[user_group.value]
    reqd_config={x:config[x] for x in config_options} 
    flattened_options=flatten_dict(reqd_config)

    updated_config=reqd_config.copy()
    for option in flattened_options:
        updated_value=typer.prompt(f"Enter new value for {option} [Press Enter to skip]",default=flattened_options[option])
        updated_config[option]=updated_value

    config_save_path=typer.prompt(f"Enter name for updated config file [Press Enter to skip]", default="config.cli_updated.yaml")
    
    save_config_file(config_save_path,unflatten_dict(updated_config))

@app.command()
def main():
    # Display the welcome message
    typer.echo(f"[bold green]🚀 Welcome to PyPSA-Earth config CLI application!")

if __name__ == "__main__":

    
    app()