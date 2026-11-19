# SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Tests for sector-resolved district-heating demand."""

import unittest

import numpy as np
import pandas as pd

from scripts.district_heat_loads import split_urban_central_heat_demand


class TestUrbanCentralHeatLoads(unittest.TestCase):
    def setUp(self):
        self.nodes = pd.Index(["JP0", "JP1"])
        self.snapshots = pd.date_range("2023-01-01", periods=3, freq="h")

        sectors = [
            "residential water",
            "residential space",
            "services water",
            "services space",
        ]

        columns = pd.MultiIndex.from_product([sectors, self.nodes])

        self.heat_demand = pd.DataFrame(
            [
                [2, 3, 10, 12, 1, 2, 4, 5],
                [3, 4, 15, 18, 2, 3, 5, 6],
                [4, 5, 12, 14, 3, 4, 6, 7],
            ],
            index=self.snapshots,
            columns=columns,
        )

        self.district_fraction = pd.Series(
            [0.3, 0.6],
            index=self.nodes,
        )

    def test_preserves_original_hourly_demand(self):
        loss = 0.15

        loads = split_urban_central_heat_demand(
            self.heat_demand,
            self.nodes,
            self.district_fraction,
            loss,
        )

        original = (
            self.heat_demand.T.groupby(level=1)
            .sum()
            .T[self.nodes]
            .multiply(self.district_fraction * (1 + loss))
        )

        reconstructed = sum(loads.values())

        np.testing.assert_allclose(
            reconstructed.to_numpy(),
            original.to_numpy(),
            rtol=1e-12,
            atol=1e-10,
        )

        self.assertEqual(
            set(loads),
            {"residential", "services", "distribution losses"},
        )

        self.assertTrue((loads["distribution losses"] >= 0).all().all())

    def test_zero_distribution_losses(self):
        loads = split_urban_central_heat_demand(
            self.heat_demand,
            self.nodes,
            self.district_fraction,
            0.0,
        )

        self.assertTrue((loads["distribution losses"] == 0).all().all())


if __name__ == "__main__":
    unittest.main()
