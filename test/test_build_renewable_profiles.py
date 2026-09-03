# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import unittest
from types import SimpleNamespace

import build_renewable_profiles


class TestCheckCutoutMatch(unittest.TestCase):
    def setUp(self):
        self.cutout = SimpleNamespace(bounds=(0, 0, 10, 10))

    def test_no_overlap_reports_both_bounds(self):
        regions = SimpleNamespace(total_bounds=(20, 2, 30, 8))

        with self.assertRaises(AssertionError) as raised:
            build_renewable_profiles.check_cutout_match(self.cutout, regions)

        message = str(raised.exception)
        self.assertIn("The requested region does not overlap the cutout.", message)
        self.assertIn("Cutout bounds: (0, 0, 10, 10)", message)
        self.assertIn("Requested region bounds: (20, 2, 30, 8)", message)
        self.assertIn("exceeds cutout bounds on: east", message)

    def test_partial_coverage_reports_uncovered_sides(self):
        regions = SimpleNamespace(total_bounds=(-2, 2, 8, 12))

        with self.assertLogs(
            build_renewable_profiles.logger, level="WARNING"
        ) as captured:
            build_renewable_profiles.check_cutout_match(self.cutout, regions)

        message = "\n".join(captured.output)
        self.assertIn(
            "Weather data only partially covers the requested region.", message
        )
        self.assertIn("Cutout bounds: (0, 0, 10, 10)", message)
        self.assertIn("Requested region bounds: (-2, 2, 8, 12)", message)
        self.assertIn("exceeds cutout bounds on: west, north", message)


if __name__ == "__main__":
    unittest.main()
