# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import os
import yaml


GEGIS_YEARS = [2013, 2015, 2018]
DEMCAST_YEARS = list(range(1999, 2025))


def check_load_weather_year(config):
    issues = []

    source = config.get("load_options", {}).get("source")
    weather_year = config.get("load_options", {}).get("weather_year")

    if source == "gegis" and weather_year not in GEGIS_YEARS:
        issues.append(
            {
                "level": "error",
                "check": "load_options.weather_year",
                "message": (
                    f"Requested load_options.weather_year={weather_year} is not "
                    f"available for load_options.source='gegis'. Available weather "
                    f"years are: {GEGIS_YEARS}."
                ),
            }
        )

    if source == "demcast" and weather_year not in DEMCAST_YEARS:
        issues.append(
            {
                "level": "error",
                "check": "load_options.weather_year",
                "message": (
                    f"Requested load_options.weather_year={weather_year} is not "
                    f"available for load_options.source='demcast'. Available weather "
                    f"years are: {DEMCAST_YEARS}."
                ),
            }
        )

    if source not in ["gegis", "demcast"]:
        issues.append(
            {
                "level": "error",
                "check": "load_options.source",
                "message": (
                    f"Unknown load_options.source='{source}'. Expected one of: "
                    "'gegis', 'demcast'."
                ),
            }
        )

    return issues


def check_config_consistency(config):
    issues = []

    issues.extend(check_load_weather_year(config))

    return issues


def write_config_consistency_report(config, output_path):
    issues = check_config_consistency(config)

    report = {
        "status": "ok" if not issues else "issues_found",
        "issues": issues,
    }

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w") as f:
        yaml.safe_dump(report, f, sort_keys=False)

    return issues