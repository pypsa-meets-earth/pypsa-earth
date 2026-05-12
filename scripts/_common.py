# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

from functools import lru_cache
from pathlib import Path

import pandas as pd


@lru_cache
def load_data_versions(file_path) -> pd.DataFrame:
    """Load the dataset versions from a CSV file and process the tags into individual columns.

    Arguments
    ---------
    file_path : str
        The path to the CSV file containing dataset version information.
    Returns
    -------
    pd.DataFrame
        A DataFrame containing the dataset version information with tags as individual boolean columns.
    """
    data_versions = pd.read_csv(
        file_path,
        dtype=str,
        na_filter=False,
        delimiter=",",
        comment="#",
    )

    # Turn space-separated tags into individual columns
    data_versions["tags"] = data_versions["tags"].str.split()
    rows_exploded_by_tags = data_versions.explode("tags")
    tags_dummies = pd.get_dummies(rows_exploded_by_tags["tags"], dtype=bool)
    enabled_tags_matrix = tags_dummies.groupby(tags_dummies.index).max()
    data_versions = data_versions.join(enabled_tags_matrix)

    return data_versions


def dataset_version(
    name: str, config: dict, versions_file_path: str = "data/versions.csv"
) -> pd.Series:
    """
    Return the dataset version information and url for a given dataset name.

    The dataset name is used to determine the source and version of the dataset from the configuration.
    Then the 'data/versions.csv' file is queried to find the matching dataset entry.

    Parameters:
    name: str
        The name of the dataset to retrieve version information for.
    config: dict
        The configuration dictionary containing dataset information, including source and version.
    versions_file_path: str, optional
        The path to the CSV file containing dataset version information (default is 'data/versions.csv').

    Returns:
    pd.Series
        A pandas Series containing the dataset version information, including source, version, tags, and URL
    """

    dataset_config = config["data"][name]

    # To use PyPSA-Zambia as a snakemake module, the path to the versions.csv file needs to be
    # registered relative to the current file with Snakemake:

    data_versions = load_data_versions(versions_file_path)

    dataset = data_versions.loc[
        (data_versions["dataset"] == name)
        & (data_versions["source"] == dataset_config["source"])
        & (data_versions["supported"])  # Limit to supported versions only
        & (
            data_versions["version"] == dataset_config["version"]
            if "latest" != dataset_config["version"]
            else True
        )
        & (data_versions["latest"] if "latest" == dataset_config["version"] else True)
    ]

    if dataset.empty:
        raise ValueError(
            f"Dataset '{name}' with source '{dataset_config['source']}' for '{dataset_config['version']}' not found in data/versions.csv."
        )

    # Return single-row DataFrame as a Series
    dataset = dataset.squeeze()

    # Generate output folder path in the `data` directory
    dataset["folder"] = Path(
        "data", name, dataset["source"], dataset["version"]
    ).as_posix()

    return dataset
