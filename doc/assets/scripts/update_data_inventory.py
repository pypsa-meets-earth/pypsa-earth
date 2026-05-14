#!/usr/bin/env python3
"""
Extract metadata for the datasets used from `data/versions.csv` and
add them into `doc/data_inventory.csv`

How it runs:
    - Automatically via the 'update-data-inventory' pre-commit hook
    - Manually: python doc/assets/scripts/extract_config_snippets.py
"""

import csv
from pathlib import Path

import pandas as pd


def update_inventory(versions_csv="versions.csv", inventory_csv="data_inventory.csv"):
    """
    Extract meta-data on the latest datasets available to run the workflow
    from `versions_csv` and place them into `inventory_csv`
    """
    priority_order = {"primary": 0, "build": 1, "archive": 2}
    category_type = pd.CategoricalDtype(categories=priority_order, ordered=True)

    versions_df = pd.read_csv(versions_csv)

    versions_df = (
        versions_df[versions_df["source"].str.lower().isin(priority_order.keys())]
        .assign(priority=lambda d: d["source"].str.lower().astype(category_type))
        .sort_values(["dataset", "priority"])
        .drop_duplicates("dataset", keep="first")[["dataset", "url"]]
    )

    inventory_df = pd.read_csv(inventory_csv)
    missed_versions_df = versions_df.loc[
        ~versions_df["dataset"]
        .str.lower()
        .isin(inventory_df["Short name"].str.lower()),
        ["dataset", "url"],
    ]

    if missed_versions_df.empty:
        return

    rows_to_append = missed_versions_df.rename(
        columns={"dataset": "Short name", "url": "Link to website"}
    ).assign(**{"Long name": "", "Description": "", "Owner": "", "License": ""})[
        [
            "Short name",
            "Long name",
            "Description",
            "Owner",
            "Link to website",
            "License",
        ]
    ]

    rows_to_append.to_csv(
        inventory_csv,
        mode="a",
        index=False,
        header=False,
        quoting=csv.QUOTE_MINIMAL,
    )


def main():
    versions_file = Path("data", "versions.csv")
    inventory_file = Path("doc", "data_inventory.csv")

    update_inventory(versions_csv=versions_file, inventory_csv=inventory_file)


if __name__ == "__main__":
    main()
