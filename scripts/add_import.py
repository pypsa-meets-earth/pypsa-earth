import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa

logger = logging.getLogger(__name__)
from _helpers import configure_logging


def load_h2_import_nodes(path):
    """
    Load prepared H2 import nodes.

    Required columns
    ----------------
    bus : str
        Base bus name without the H2 suffix.
    p_nom_max_pipeline : float
        Maximum pipeline-based H2 import capacity at the node.
    p_nom_max_port : float
        Maximum port-based H2 import capacity at the node.

    Optional columns
    ----------------
    p_nom_max_total

    Parameters
    ----------
    path : str or path-like
        Path to the prepared H2 import nodes CSV.

    Returns
    -------
    pandas.DataFrame
        Cleaned dataframe with the required columns.
    """
    df = pd.read_csv(
        path,
        index_col=None,
        keep_default_na=False,
    ).copy()

    required_columns = {"bus", "p_nom_max_pipeline", "p_nom_max_port"}
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(
            f"Missing required columns in H2 import nodes file {path}: "
            f"{sorted(missing_columns)}"
        )

    df["bus"] = df["bus"].astype(str).str.strip()
    df["p_nom_max_pipeline"] = pd.to_numeric(df["p_nom_max_pipeline"], errors="coerce")
    df["p_nom_max_port"] = pd.to_numeric(df["p_nom_max_port"], errors="coerce")

    df = df.dropna(subset=["bus", "p_nom_max_pipeline", "p_nom_max_port"])
    df = df[df["bus"] != ""]

    if (df["p_nom_max_pipeline"] < 0).any():
        invalid = df.loc[df["p_nom_max_pipeline"] < 0, ["bus", "p_nom_max_pipeline"]]
        raise ValueError(
            "Negative values found in 'p_nom_max_pipeline':\n"
            f"{invalid.to_string(index=False)}"
        )

    if (df["p_nom_max_port"] < 0).any():
        invalid = df.loc[df["p_nom_max_port"] < 0, ["bus", "p_nom_max_port"]]
        raise ValueError(
            "Negative values found in 'p_nom_max_port':\n"
            f"{invalid.to_string(index=False)}"
        )

    return df


def get_h2_target_buses(n, import_nodes):
    """
    Map prepared import nodes to existing hydrogen buses in the network.

    Parameters
    ----------
    n : pypsa.Network
        Network with H2 buses already present.
    import_nodes : pandas.DataFrame
        Prepared import nodes dataframe with base bus names.

    Returns
    -------
    pandas.DataFrame
        Copy of dataframe with an added column ``h2_bus``.
        Rows without an existing H2 bus are removed.
    """
    df = import_nodes.copy()

    df["h2_bus"] = df["bus"].astype(str) + " H2"
    df = df[df["h2_bus"].isin(n.buses.index)].copy()

    return df


def add_pipeline_h2_imports(n, import_nodes):
    """
    Add pipeline-based H2 import generators.

    Pipeline imports are modeled as extendable generators with
    ``p_nom_max = p_nom_max_pipeline``.
    """
    sources = snakemake.params.imports_config["sources"]
    if "pipelines" not in sources:
        return

    df = import_nodes.copy()
    df = df[df["p_nom_max_pipeline"] > 0].copy()

    if df.empty:
        return

    marginal_cost = snakemake.params.imports_config["price"]["H2_pipe"]

    names = df["bus"].astype(str) + " import H2 pipeline"

    n.madd(
        "Generator",
        names,
        bus=df["h2_bus"],
        carrier="import H2 pipeline",
        p_nom_extendable=True,
        p_nom_max=df["p_nom_max_pipeline"],
        marginal_cost=marginal_cost,
    )


def add_port_h2_imports(n, import_nodes):
    """
    Add port-based H2 import generators.

    Port imports are modeled as extendable generators with
    ``p_nom_max = p_nom_max_port``.
    """
    sources = snakemake.params.imports_config["sources"]
    if "ports" not in sources:
        return

    df = import_nodes.copy()
    df = df[df["p_nom_max_port"] > 0].copy()

    if df.empty:
        return

    marginal_cost = snakemake.params.imports_config["price"]["H2_ship"]

    names = df["bus"].astype(str) + " import H2 port"

    n.madd(
        "Generator",
        names,
        bus=df["h2_bus"],
        carrier="import H2 port",
        p_nom_extendable=True,
        p_nom_max=df["p_nom_max_port"],
        marginal_cost=marginal_cost,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "add_import",
            simpl="",
            clusters="10",
        )

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)

    import_nodes = load_h2_import_nodes(snakemake.input.h2_import_nodes)

    import_nodes = get_h2_target_buses(n,import_nodes,)

    add_pipeline_h2_imports(n,import_nodes,)

    add_port_h2_imports(n,import_nodes,)

    n.export_to_netcdf(snakemake.output[0])
    
    logger.info("Network successfully exported")
    