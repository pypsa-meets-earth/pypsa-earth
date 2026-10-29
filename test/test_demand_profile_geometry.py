# SPDX-FileCopyrightText: 2026 PyPSA-Earth Contributors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from build_demand_profiles import shapes_to_shapes


def test_shapes_to_shapes_with_nonconsecutive_indices():
    orig_geoms = [
        box(0, 0, 1, 1),
        box(1, 0, 2, 1),
    ]
    dest_geoms = [
        box(0, 0, 2, 1),
        box(0, 0, 1, 1),
    ]

    orig = gpd.GeoSeries(orig_geoms, index=[10, 20])
    dest = gpd.GeoSeries(dest_geoms, index=[30, 40])

    expected = np.array(
        [
            [0.5, 0.5],
            [1.0, 0.0],
        ]
    )

    result = shapes_to_shapes(orig, dest).toarray()

    np.testing.assert_allclose(result, expected)


def test_shapes_to_shapes_with_lists():
    orig = [
        box(0, 0, 1, 1),
        box(1, 0, 2, 1),
    ]
    dest = [
        box(0, 0, 2, 1),
        box(0, 0, 1, 1),
    ]

    expected = np.array(
        [
            [0.5, 0.5],
            [1.0, 0.0],
        ]
    )

    result = shapes_to_shapes(orig, dest).toarray()

    np.testing.assert_allclose(result, expected)
