# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Command-line interface (CLI) for inspecting and retrieving databundles.

This non-workflow helper utility scans configured databundles, checks whether
the expected local files exist on disk, presents an interactive checklist using
the ``rich`` library, and provides options to retrieve missing bundles either
automatically or manually via cloud storage links (Zenodo, Google Drive, direct).

Typical usage
-------------
Run the script directly from the repository root:

```bash
python scripts/non_workflow/databundle_cli.py
```

Relevant Settings
-----------------
The script reads configuration settings from ``config.yaml`` / ``config.default.yaml``
and ``configs/bundle_config.yaml``:

```yaml
enable:
  progress_bar: true

databundles:
  # Bundle definitions with URLs and output paths
```

Inputs
------
- ``configs/bundle_config.yaml``: Databundle definitions, source URLs, and expected output paths.
- ``logs/databundle_cli.yaml`` (optional): Custom parameters logged by previous workflow runs.

Outputs
-------
- Downloaded and extracted dataset files into local project directories (e.g. ``data/``, ``cutouts/``).
"""

import os
import sys
import textwrap
from typing import Any

import pandas as pd
import yaml

# Add parent directory to sys.path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, parent_dir)

from _helpers import mock_snakemake
from retrieve_databundle_light import (
    load_databundle_config,
    retrieve_databundle,
)
from rich.console import Console
from rich.markdown import Markdown
from rich.table import Table


def console_markdown(markdown: str, lvl: int = 1) -> None:
    """
    Print markdown text to the console with optional indentation.

    Parameters
    ----------
    markdown : str
        The markdown text to be printed.
    lvl : int, default 1
        The level of indentation (number of times to dedent the markdown text).

    Returns
    -------
    None
    """
    console = Console()
    for _ in range(lvl):
        markdown = textwrap.dedent(markdown)
    md = Markdown(markdown)

    console.print(md)


def console_table(
    dataframe: pd.DataFrame, table_kw: dict[str, Any] | None = None
) -> None:
    """
    Print a DataFrame as a rich table in the console.

    Parameters
    ----------
    dataframe : pd.DataFrame
        The DataFrame to be printed as a table.
    table_kw : dict[str, Any], optional
        Keyword arguments for customizing the table appearance.

    Returns
    -------
    None
    """
    if table_kw is None:
        table_kw = {}

    df = dataframe.reset_index()

    # Initialize the console
    console = Console()

    # Create a rich table
    table = Table(show_header=True, **table_kw)

    # Add columns from the DataFrame
    for column in df.columns:
        table.add_column(str(column), overflow="fold")

    # Add rows from the DataFrame
    for row in df.itertuples(index=False):
        if hasattr(row, "retrieve"):
            color = "green" if getattr(row, "retrieve", False) else "red"
            table.add_row(*map(str, row), style=color)
        else:
            table.add_row(*map(str, row))

    console.print(table)


def databundle_check(
    bundles_to_download: list[str], config_databundles: dict[str, Any]
) -> list[str]:
    """
    Check the status of the files in the databundle and print a checklist to the console.

    Parameters
    ----------
    bundles_to_download : list[str]
        A list of databundle names to check.
    config_databundles : dict[str, Any]
        A dictionary containing the configuration for each databundle.

    Returns
    -------
    list[str]
        A list of missing databundle names.
    """
    df_file = pd.DataFrame()

    for bundlename in bundles_to_download:
        df = (
            pd.DataFrame.from_dict(config_databundles[bundlename]["output"])
            .rename(columns={0: "filepath"})
            # .assign(category=config_databundles[bundlename]["category"])
            .assign(databundle=bundlename)
            .set_index("filepath")
        )

        # Filter out rows where the 'file path' contains '*'
        df = df[~df.index.str.contains(r"\*", na=False)]
        df = df[~df.index.str.endswith("/")]

        df_file = pd.concat([df_file, df])

    df_file["retrieve"] = [os.path.isfile(str(i)) for i in df_file.index]

    console_table(
        df_file, table_kw={"title": "Databundle Checklist"}
    )  # [["databundle","retrieve"]]

    missing_bundles = list(df_file.loc[~df_file.retrieve, "databundle"].unique())

    if len(missing_bundles) == 0:
        console_markdown("No data is missing")

    else:
        extractable_bundle = [
            bundlename
            for bundlename in missing_bundles
            if set(config_databundles[bundlename]["urls"].keys()).intersection(
                ["zenodo", "gdrive", "direct"]
            )
        ]

        console_markdown(
            f"""It looks like there are missing files from: **{", ".join(missing_bundles)}**"""
        )

        if extractable_bundle:
            df_link = pd.DataFrame()

            for bundlename in extractable_bundle:
                df = (
                    pd.DataFrame(config_databundles[bundlename]["urls"].items())
                    .rename(columns={0: "host", 1: "urls"})
                    .assign(databundle=bundlename)
                )

                df.loc[df.index[1:], "databundle"] = ""
                df = df.set_index("databundle")

                df_link = pd.concat([df_link, df])

            console_table(df_link, table_kw={"title": "Retrievable Databundle"})

            console_markdown(
                f"""The databundle **{", ".join(extractable_bundle)}** can be retrieved either automatically or manually downloaded."""
            )

        unextractable_bundle = [
            bundlename
            for bundlename in missing_bundles
            if bundlename not in extractable_bundle
        ]

        if unextractable_bundle:
            console_markdown(
                f"""The databundle **{", ".join(unextractable_bundle)}** can only be retrieved automatically."""
            )

    console_markdown(f"""
    Options:

    - **check**: update the checklist table to see if the file is now included
    - **all**: get all missing databundles: **{", ".join(missing_bundles)}**
    - **rerun**: get all databundles again: **{", ".join(bundles_to_download)}**
    - **bundle_...**: get selected databundles, can be more than one
    - Press **ENTER** to exit
    """)

    return missing_bundles


if __name__ == "__main__":

    console_markdown("""

    # PyPSA-Earth Databundle Retrieval Command-Line Interface (CLI)

    Find the missing file for your snakemake run and solve them your own way.
    To prevent snakemake from trying to re-download the file again during the workflow, make sure to:

    - Set `retrieve_databundle` as false

    """)

    log_path = "logs/databundle_cli.yaml"
    rootpath = "."

    # Load the rule from YAML if it exists, otherwise use default
    if os.path.exists(log_path):
        with open(log_path, encoding="utf-8") as f:
            snakemake_rule = yaml.safe_load(f)
    else:
        snakemake_rule = {"rulename": "retrieve_databundle_light"}

    snakemake = mock_snakemake(**snakemake_rule)

    if "params" in snakemake_rule:
        snakemake.params = snakemake_rule["params"]

    # load enable configuration
    config_enable = snakemake.config["enable"]
    # load databundle configuration
    config_bundles = load_databundle_config(snakemake.config["databundles"])
    disable_progress = not config_enable["progress_bar"]
    hydrobasins_level = snakemake.params["hydrobasins_level"]

    bundles_to_download = snakemake.params["bundles_to_download"]
    config_databundles = snakemake.config["databundles"]

    while True:
        missing_bundles = databundle_check(bundles_to_download, config_databundles)
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
            console_markdown(f"**{', '.join(listed_bundles)}** selected")
            bundles = tuple(listed_bundles)
        else:
            console_markdown("incorrect input")
            continue

        retrieve_databundle(
            bundles,
            config_bundles,
            hydrobasins_level,
            rootpath=rootpath,
            disable_progress=disable_progress,
        )

    if os.path.exists(log_path):
        os.remove(log_path)
        print(f"{log_path} deleted.")
