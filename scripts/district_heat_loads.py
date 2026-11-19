# SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Build sector-resolved district-heating demand profiles."""

import pandas as pd


def split_urban_central_heat_demand(
    heat_demand: pd.DataFrame,
    nodes: pd.Index,
    district_fraction: pd.Series,
    distribution_loss: float,
) -> dict[str, pd.DataFrame]:
    """Preserve sectoral final demand and represent distribution losses separately."""
    loads = {}

    for sector in ("residential", "services"):
        sector_demand = (
            heat_demand[[f"{sector} water", f"{sector} space"]]
            .T.groupby(level=1)
            .sum()
            .T[nodes]
        )

        loads[sector] = sector_demand.multiply(district_fraction.loc[nodes])

    total_delivered = loads["residential"] + loads["services"]

    loads["distribution losses"] = total_delivered * distribution_loss

    return loads
