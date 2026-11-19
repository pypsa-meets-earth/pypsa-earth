# SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Utilities for selecting sector technology cost years."""

import warnings
from collections.abc import Mapping
from pathlib import Path

import yaml


def read_available_cost_years(metadata_path: str | Path) -> list[int]:
    """Read available cost years from technology-data metadata."""
    with open(metadata_path, encoding="utf-8") as f:
        metadata = yaml.safe_load(f)

    years = metadata.get("years") if isinstance(metadata, dict) else None

    if not years:
        raise ValueError(
            "No cost years found in the configured technology-data release."
        )

    return sorted({int(year) for year in years})


def resolve_sector_cost_year(
    planning_horizon: int,
    *,
    foresight: str,
    costs: Mapping,
    retrieve_cost_data: bool,
    available_years: list[int] | None = None,
) -> int:
    """Resolve the sector cost year according to the configured policy."""
    if foresight == "overnight":
        return int(costs["year"])

    horizon = int(planning_horizon)
    selection = costs.get("year_selection", "closest")

    if isinstance(selection, Mapping):
        mapping = {int(year): int(cost_year) for year, cost_year in selection.items()}

        if horizon not in mapping:
            raise ValueError(
                f"No cost year defined for planning horizon {horizon} "
                "in costs.year_selection."
            )

        return mapping[horizon]

    if selection != "closest":
        raise ValueError(
            "costs.year_selection must be 'closest' or an explicit mapping."
        )

    if not retrieve_cost_data:
        return horizon

    if available_years is None:
        raise ValueError(
            "Available cost years are required for closest-year selection."
        )

    years = sorted({int(year) for year in available_years})

    if not years:
        raise ValueError("No available technology-data cost years provided.")

    selected_year = min(
        years,
        key=lambda year: (abs(year - horizon), year),
    )

    if selected_year != horizon:
        warnings.warn(
            f"Cost year {horizon} is not available in technology-data. "
            f"Using closest available year {selected_year} instead.",
            stacklevel=2,
        )

    return selected_year
