import textwrap
from rich.table import Table
from rich.console import Console
from rich.markdown import Markdown
from retrieve_databundle_light import *
from _helpers import create_logger

logger = create_logger(__name__)

import yaml

configfile = ["config.default.yaml",
              "configs/bundle_config.yaml",
              "config.yaml"]

config = {}
for c in configfile:
    if os.path.isfile(c):
        with open(c) as file:
            config_append = yaml.safe_load(file)

        config.update(config_append)

def console_markdown(markdown, lvl=1):

    console = Console()
    for i in range(lvl):
        markdown = textwrap.dedent(markdown)
    md = Markdown(markdown)

    return console.print(md)

def console_table(dataframe, table_kw={}):
    df = dataframe.reset_index()

    # Initialize the console
    console = Console()
    
    # Create a rich table
    table = Table(show_header=True, **table_kw)
    
    # Add columns from the DataFrame
    for column in df.columns:
        table.add_column(column, overflow="fold")
    
    # Add rows from the DataFrame
    for row in df.itertuples(index=False):
        if hasattr(row, 'retrieve'):
            color = "green" if getattr(row, "retrieve", False) else "red"
            table.add_row(*map(str, row), style=color)
        else:
            table.add_row(*map(str, row))
        
    # Print the table to the console
    return console.print(table)


def databundle_check(bundles_to_download):

    df_file = pd.DataFrame()

    for bundlename in bundles_to_download:
        df = (pd.DataFrame.from_dict(config["databundles"][bundlename]["output"])
            .rename(columns={0:"filepath"})
            #.assign(category=config["databundles"][bundlename]["category"])
            .assign(databundle=bundlename)
            .set_index("filepath")
            )

        # Filter out rows where the 'file path' contains '*'
        df = df[~df.index.str.contains(r'\*', na=False)]
        df = df[~df.index.str.endswith("/")]
        
        df_file = pd.concat([df_file,df])

    df_file["retrieve"] = [os.path.isfile(f"{i}") for i in df_file.index]

    console_table(df_file, table_kw={"title":"Databundle Checklist"}) #[["databundle","retrieve"]]

    missing_bundles = list(df_file.loc[~df_file.retrieve, "databundle"].unique())

    if len(missing_bundles) == 0:
        console_markdown("No data is missing")
    
    else:
        extractable_bundle = [
            bundlename for bundlename in missing_bundles
            if set(config["databundles"][bundlename]["urls"].keys()).intersection(["zenodo", "gdrive", "direct"])
        ]

        console_markdown(f"""It looks like there are missing files from: **{", ".join(missing_bundles)}**""")

        if extractable_bundle:
            df_link = pd.DataFrame()

            for bundlename in extractable_bundle:
                df = (pd.DataFrame(config["databundles"][bundlename]["urls"].items())
                    .rename(columns={0:"host", 1: "urls"})
                    .assign(databundle=bundlename)
                    )
                
                df.loc[df.index[1:], "databundle"] = ""
                df = df.set_index("databundle")
                
                df_link = pd.concat([df_link,df])

            console_table(df_link, table_kw={"title":"Retrievable Databundle"})
            
            console_markdown(f"""The databundle **{", ".join(extractable_bundle)}** can be retrieved either automatically or manually downloaded.""")

        unextractable_bundle = [bundlename for bundlename in missing_bundles if bundlename not in extractable_bundle]

        if unextractable_bundle:
            console_markdown(f"""The databundle **{", ".join(unextractable_bundle)}** can only be retrieved automatically.""")

    console_markdown(f"""
    Options:

    - **check**: update the Databundle Checklist to see if the file is included
    - **all**: retrieve all missing databundles, namely **{", ".join(missing_bundles)}**
    - **rerun**: retrieve all databundles again, namely **{", ".join(bundles_to_download)}**
    - **bundle_...**: retrieve the selected databundles, can be more than one
    - Press **ENTER** to end this loop
    """)

    return missing_bundles

if __name__ == "__main__":

    console_markdown(""" 

    # PyPSA-Earth Databundle Retrival Command-Line Interface (CLI)

    Find the missing file for your snakemake run and solve them your own way. Please make sure to:

    - Set `retrieve_databundle` as false

    """)

    rootpath = "."
    tutorial = config["tutorial"]
    countries = config["countries"]
    config_enable = config["enable"]
    disable_progress = not config_enable["progress_bar"]
    hydrobasins_level = config["renewable"]["hydro"]["hydrobasins_level"]

    config_bundles = load_databundle_config(config["databundles"])

    bundles_to_download = get_best_bundles(
            countries, config_bundles, tutorial, config_enable
        )    
    
    while True:
        missing_bundles = databundle_check(bundles_to_download)
        answer = input("Input: ")

        if not answer:
            break

        answer = answer.split(" ")
        listed_bundles = [item for item in answer if item in bundles_to_download]
        if "check" in answer:
            console_markdown("**check** is selected")
            continue
        elif "all" in answer:
            console_markdown("**all** is selected")
            bundles = missing_bundles
        elif "rerun" in answer:
            console_markdown("**rerun** is selected")
            bundles = bundles_to_download
        elif listed_bundles:
            console_markdown(f"**{listed_bundles[0]}** is selected")
            bundles = listed_bundles[0], 
        else:
            console_markdown("incorrect input")
            continue
            
        retrieve_databundle(bundles, 
                config_bundles,
                hydrobasins_level, 
                rootpath=rootpath, 
                disable_progress=disable_progress
                )