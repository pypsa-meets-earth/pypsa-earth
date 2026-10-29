# SPDX-FileCopyrightText: 2026 PyPSA-Earth Contributors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import pandas as pd


def test_hydropower_inflow_aggregation():
    snapshots = pd.date_range("2030-01-01", periods=3, freq="h")

    inflow_agg = pd.DataFrame(
        {"bus_A hydro": [10.0, 20.0, 30.0]},
        index=snapshots,
    )

    inflow_agg_ror = pd.DataFrame(
        {"bus_B ror": [2.0, 3.0, 4.0]},
        index=snapshots,
    )

    inflow_agg_hydro = pd.DataFrame(
        {
            "bus_A hydro": [1.0, 2.0, 3.0],
            "bus_C hydro": [5.0, 6.0, 7.0],
        },
        index=snapshots,
    )

    result = (
        pd.concat([inflow_agg, inflow_agg_ror, inflow_agg_hydro], axis=1)
        .T.groupby(level=0)
        .sum()
        .T
    )

    expected = pd.DataFrame(
        {
            "bus_A hydro": [11.0, 22.0, 33.0],
            "bus_B ror": [2.0, 3.0, 4.0],
            "bus_C hydro": [5.0, 6.0, 7.0],
        },
        index=snapshots,
    )

    pd.testing.assert_frame_equal(
        result.sort_index(axis=1),
        expected.sort_index(axis=1),
    )
