# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Check that the shapes built in build_shapes keep their CRS.

GADM and EEZ inputs are replaced by small synthetic geometries, so these
tests do not download any data.
"""

import sys
import warnings

import geopandas as gpd
import pytest
from shapely.geometry import MultiPolygon, box

sys.path.append("./scripts")

import build_shapes  # noqa: E402

GEO_CRS = "EPSG:4326"


@pytest.fixture
def fake_gadm(monkeypatch):
    def _get_GADM_layer(countries, layer_id, geo_crs, *args, **kwargs):
        return gpd.GeoDataFrame(
            {"GID_0": ["AA", "BB"]},
            geometry=[
                MultiPolygon([box(0, 0, 2, 2), box(5, 5, 5.5, 5.5)]),
                box(2, 0, 4, 2),
            ],
            crs=geo_crs,
        )

    monkeypatch.setattr(build_shapes, "get_GADM_layer", _get_GADM_layer)


@pytest.fixture
def fake_eez(monkeypatch):
    def _load_EEZ(countries, geo_crs, EEZ_gpkg):
        return gpd.GeoDataFrame(
            {"name": ["AA", "BB"]},
            geometry=[box(-1, -1, 1, 3), box(3, -1, 5, 3)],
            crs=geo_crs,
        )

    monkeypatch.setattr(build_shapes, "load_EEZ", _load_EEZ)


def _crs_warnings(record):
    return [str(w.message) for w in record if "CRS" in str(w.message)]


def test_countries_keeps_crs(fake_gadm):
    country_shapes = build_shapes.countries(["AA", "BB"], GEO_CRS, "set_by_country")

    assert country_shapes.crs == GEO_CRS


def test_eez_keeps_crs_without_warnings(fake_gadm, fake_eez):
    country_shapes = build_shapes.countries(["AA", "BB"], GEO_CRS, "set_by_country")

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        offshore_shapes = build_shapes.eez(
            ["AA", "BB"], GEO_CRS, country_shapes, "unused.gpkg"
        )

    assert offshore_shapes.crs == GEO_CRS
    assert _crs_warnings(record) == []
    assert not offshore_shapes.empty


def test_eez_without_simplification_keeps_crs(fake_gadm, fake_eez):
    country_shapes = build_shapes.countries(["AA", "BB"], GEO_CRS, "set_by_country")

    offshore_shapes = build_shapes.eez(
        ["AA", "BB"], GEO_CRS, country_shapes, "unused.gpkg", simplify_gadm=False
    )

    assert offshore_shapes.crs == GEO_CRS


def test_eez_buffer_distance_still_removes_coastal_strip(fake_gadm, fake_eez):
    country_shapes = build_shapes.countries(["AA", "BB"], GEO_CRS, "set_by_country")

    with warnings.catch_warnings():
        # buffering in degrees is expected to warn about the geographic CRS
        warnings.simplefilter("ignore", UserWarning)
        no_buffer = build_shapes.eez(
            ["AA", "BB"], GEO_CRS, country_shapes, "unused.gpkg"
        )
        with_buffer = build_shapes.eez(
            ["AA", "BB"], GEO_CRS, country_shapes, "unused.gpkg", distance=0.1
        )

    assert with_buffer.crs == GEO_CRS
    assert (
        with_buffer.to_crs("EPSG:3857").area < no_buffer.to_crs("EPSG:3857").area
    ).all()
