# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: PyPSA-ASEAN, PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Add transmission projects to the network.
"""

from pathlib import Path

import pandas as pd
import pypsa
from _helpers import configure_logging, create_logger

logger = create_logger(__name__)


def attach_transmission_projects(
    n: pypsa.Network, transmission_projects: list[str]
) -> None:
    logger.info("Adding transmission projects to network.")
    for path in transmission_projects:
        path = Path(path)
        logger.info(path)
        df = pd.read_csv(path, index_col=0, dtype={"bus0": str, "bus1": str})
        if df.empty:
            continue

        if "new_buses" in path.name:
            logger.info(f"Adding new buses from {path}")

            missing_coords = df[df["x"].isnull() | df["y"].isnull()]
            if not missing_coords.empty:
                logger.warning(
                    f"Dropping {len(missing_coords)} buses from {path} "
                    f"because they have missing coordinates:\n"
                    f"{missing_coords.index.tolist()}"
                )
                # Drop rows with missing coordinates
                df = df.drop(missing_coords.index)

            if not df.empty:
                df.rename(columns={"x": "lon", "y": "lat"}, inplace=True)
                n.madd("Bus", df.index, **df.to_dict(orient="list"))
            else:
                logger.warning(f"No buses added from {path} after filtering.")

        elif "new_lines" in path.name:
            logger.info(f"Adding new lines from {path}")
            n.madd("Line", df.index, **df.to_dict(orient="list"))

        elif "new_links" in path.name:
            logger.info(f"Adding new links from {path}")
            n.madd("Link", df.index, **df.to_dict(orient="list"))

        elif "adjust_lines" in path.name:
            logger.info(f"Adjusting lines from {path}")
            n.lines.update(df)

        elif "adjust_links" in path.name:
            logger.info(f"Adjusting links from {path}")
            n.links.update(df)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("add_transmission_projects")
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)

    params = snakemake.params

    if params["transmission_projects"]["enable"]:
        attach_transmission_projects(n, snakemake.input.transmission_projects)

    nan_buses = n.buses[n.buses[["lon", "lat"]].isna().any(axis=1)]

    if not nan_buses.empty:
        raise ValueError(
            f"ERROR: Buses with missing coordinates detected before export:\n"
            f"{nan_buses}"
        )
    else:
        logger.info("All buses have valid coordinates.")

    if {"lon", "lat"}.issubset(n.buses.columns):
        mask_xy_missing = (
            n.buses["x"].isna()
            | n.buses["y"].isna()
            | ((n.buses["x"] == 0.0) & (n.buses["y"] == 0.0))
        )

        n.buses.loc[mask_xy_missing, "x"] = n.buses.loc[mask_xy_missing, "lon"]
        n.buses.loc[mask_xy_missing, "y"] = n.buses.loc[mask_xy_missing, "lat"]

        logger.info(
            f"Patched {mask_xy_missing.sum()} buses where x/y were zero or missing "
            f"using lon/lat columns."
        )
    else:
        logger.info("No lon/lat columns in buses table; skipping patch.")

    # Get current bus index as strings
    bus_index = n.buses.index.astype(str)

    # Find buses whose names are NOT purely numeric
    non_numeric_mask = ~bus_index.str.match(r"^\d+$")
    non_numeric_buses = n.buses[non_numeric_mask]

    if not non_numeric_buses.empty:
        # Find the largest current integer bus name
        numeric_buses = n.buses[~non_numeric_mask]
        if not numeric_buses.empty:
            max_id = numeric_buses.index.astype(int).max()
        else:
            max_id = -1

        new_names = {}
        for i, old_name in enumerate(non_numeric_buses.index, start=max_id + 1):
            new_names[old_name] = str(i)

        # Rename buses in the network
        n.buses.rename(index=new_names, inplace=True)

        # Also rename any references to buses in other components
        for component in n.components:
            attr_name = component.lower() + "s"
            if hasattr(n, attr_name):
                df = getattr(n, attr_name)
                if isinstance(df, pd.DataFrame):
                    for col in df.columns:
                        if df[col].dtype == object:
                            df[col] = df[col].replace(new_names)

        logger.info(f"Renamed {len(new_names)} non-integer bus names:\n{new_names}")
    else:
        logger.info("No non-integer bus names detected.")

    # Ensure 'under_construction' column contains only boolean values
    n.links["under_construction"] = n.links["under_construction"] == True
    n.lines["under_construction"] = n.links["under_construction"] == True

    n.export_to_netcdf(snakemake.output[0])
