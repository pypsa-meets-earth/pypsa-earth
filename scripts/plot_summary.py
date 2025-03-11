# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Plots energy and cost summaries for solved networks.

Relevant Settings
-----------------
Inputs
------
Outputs
-------
Description
-----------
"""
import os

import matplotlib.pyplot as plt
import pandas as pd
from _helpers import configure_logging, create_logger

logger = create_logger(__name__)


def rename_techs(label):
    if "H2" in label:
        label = "hydrogen storage"
    elif label == "solar":
        label = "solar PV"
    elif label == "offwind-ac":
        label = "offshore wind ac"
    elif label == "offwind-dc":
        label = "offshore wind dc"
    elif label == "onwind":
        label = "onshore wind"
    elif label == "ror":
        label = "hydroelectricity"
    elif label == "hydro":
        label = "hydroelectricity"
    elif label == "PHS":
        label = "hydroelectricity"
    elif "battery" in label:
        label = "battery storage"

    return label


preferred_order = pd.Index(
    [
        "transmission lines",
        "hydroelectricity",
        "hydro reservoir",
        "run of river",
        "pumped hydro storage",
        "onshore wind",
        "offshore wind ac",
        "offshore wind dc",
        "solar PV",
        "solar thermal",
        "OCGT",
        "hydrogen storage",
        "battery storage",
    ]
)


def plot_costs(infn, snmk, fn=None):
    # For now ignore the simpl header
    cost_df = pd.read_csv(infn, index_col=list(range(3)), header=[1, 2, 3])

    df = cost_df.groupby(cost_df.index.get_level_values(2)).sum()

    # convert to billions
    df = df / 1e9

    df = df.groupby(df.index.map(rename_techs)).sum()

    to_drop = df.index[df.max(axis=1) < snmk.config["plotting"]["costs_threshold"]]

    if not to_drop.empty:
        logger.info(
            "Dropping elements from costs dataframe:\n"
            + df.loc[to_drop].to_string()
            + "\n"
        )

    df = df.drop(to_drop)

    if not to_drop.empty:
        logger.info("Remaining elements in costs dataframe:\n" + df.to_string() + "\n")

    new_index_costs = df.index.intersection(preferred_order).append(
        df.index.difference(preferred_order)
    )

    new_columns = df.sum().sort_values().index

    fig_costs, ax = plt.subplots()
    fig_costs.set_size_inches((12, 8))

    if new_index_costs.empty:
        logger.error(
            f"No costs data to plot for country {snmk.wildcards.country}, no valid plot is generated"
        )
        # create empty file to avoid issues with snakemake
        with open(fn, "w") as fp:
            pass
        return

    df.loc[new_index_costs, new_columns].T.plot(
        kind="bar",
        ax=ax,
        stacked=True,
        color=[snmk.config["plotting"]["tech_colors"][i] for i in new_index_costs],
    )

    handles, labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    ax.set_ylim([0, snmk.config["plotting"]["costs_max"]])

    ax.set_ylabel("System Cost [EUR billion per year]")

    ax.set_xlabel("")

    ax.grid(axis="y")

    ax.legend(handles, labels, ncol=4, loc="upper left")

    fig_costs.tight_layout()

    if fn is not None:
        fig_costs.savefig(fn, transparent=True)


def plot_energy(infn, snmk, fn=None):
    energy_df = pd.read_csv(infn, index_col=list(range(2)), header=[1, 2, 3])

    df = energy_df.groupby(energy_df.index.get_level_values(1)).sum()

    # convert MWh to TWh
    df = df / 1e6

    df = df.groupby(df.index.map(rename_techs)).sum()

    to_drop = df.index[
        df.abs().max(axis=1) < snmk.config["plotting"]["energy_threshold"]
    ]

    if not to_drop.empty:
        logger.info(
            "Dropping elements from energy dataframe:\n"
            + df.loc[to_drop].to_string()
            + "\n"
        )

    df = df.drop(to_drop)

    if not to_drop.empty:
        logger.info("Remaining elements in energy dataframe:\n" + df.to_string() + "\n")

    new_index_energy = df.index.intersection(preferred_order).append(
        df.index.difference(preferred_order)
    )

    new_columns = df.columns.sort_values()

    fig_energy, ax = plt.subplots()
    fig_energy.set_size_inches((12, 8))

    if new_index_energy.empty:
        logger.error(
            f"No energy data to plot for country {snmk.wildcards.country}, no valid plot is generated"
        )
        # create empty file to avoid issues with snakemake
        with open(fn, "w") as fp:
            pass
        return

    df.loc[new_index_energy, new_columns].T.plot(
        kind="bar",
        ax=ax,
        stacked=True,
        color=[snmk.config["plotting"]["tech_colors"][i] for i in new_index_energy],
    )

    handles, labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    ax.set_ylim(
        [
            snmk.config["plotting"]["energy_min"],
            snmk.config["plotting"]["energy_max"],
        ]
    )

    ax.set_ylabel("Energy [TWh/a]")

    ax.set_xlabel("")

    ax.grid(axis="y")

    ax.legend(handles, labels, ncol=4, loc="upper left")

    fig_energy.tight_layout()

    if fn is not None:
        fig_energy.savefig(fn, transparent=True)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_summary",
            summary="energy",
            network="elec",
            simpl="",
            clusters="4",
            ll="c1",
            opts="Co2L-4H",
            attr="",
            ext="pdf",
            country="all",
        )
    configure_logging(snakemake)

    summary = snakemake.wildcards.summary
    try:
        func = globals()[f"plot_{summary}"]
    except KeyError:
        logger.error(f"plotting function for {summary} has not been defined")

    func(
        os.path.join(snakemake.input[0], f"{summary}.csv"),
        snakemake,
        snakemake.output[0],
    )
