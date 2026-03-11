# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import pandas as pd
from functools import partial, lru_cache


@lru_cache
def load_data_versions(file_path):
    data_versions = pd.read_csv(
        file_path,
        dtype=str,
        na_filter=False,
        delimiter=",",
        comment="#",
    )

    # Turn space-separated tags into individual columns
    data_versions["tags"] = data_versions["tags"].str.split()
    exploded = data_versions.explode("tags")
    dummies = pd.get_dummies(exploded["tags"], dtype=bool)
    tags_matrix = dummies.groupby(dummies.index).max()
    data_versions = data_versions.join(tags_matrix)

    return data_versions


def dataset_version(
    name: str,
) -> pd.Series:
    """
    Return the dataset version information and url for a given dataset name.

    The dataset name is used to determine the source and version of the dataset from the configuration.
    Then the 'data/versions.csv' file is queried to find the matching dataset entry.

    Parameters:
    name: str
        The name of the dataset to retrieve version information for.

    Returns:
    pd.Series
        A pandas Series containing the dataset version information, including source, version, tags, and URL
    """

    dataset_config = config["data"][
        name
    ]  # TODO as is right now, it is not compatible with config_provider

    # To use PyPSA-Eur as a snakemake module, the path to the versions.csv file needs to be
    # registered relative to the current file with Snakemake:
    fp = workflow.source_path("../data/versions.csv")
    data_versions = load_data_versions(fp)

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
