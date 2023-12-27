# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Creates a dataset corresponding to the projected climate in a requested year.

Relevant Settings
-----------------

.. code:: yaml

    projection:
        climate_scenario:
        present_year:
        future_year:
        years_window:
        gcm_selection:

Inputs
------

- ``cutouts/cutout.nc``: confer :ref:`cutout`, a cutout file produced by altile
- ``data/cmip6_avr.nc``: confer :ref:`cmip6`, a global-scale dataset of the global climate
CMIP6 projections aggregated into IPCC Global Climate Atlas

Outputs
-------

- ``cutouts/proj-{cutout}.nc"``: confer :ref:`proj`, a cutout modified to account for future
climate conditions

Description
-----------

The rule :mod:`build_climate_projections` creates a cutout file which corresponds to
a requested year in the future. The parameter projections is calculated combining data
for cutout.nc which is assumed to be representative data sample for the past climate and
an ensemble of CMIP6 global climate models to account for the future climate.

A morphing approach is implemented meaning that the projection value x_future is calculated
data for a representative period present_year in the past x_past using a combination
of shift dx and stretch a:

x_future = x_past + dx + a(x_past - <x_past>|month),

where
x_future and x_past have the highest temporal resolution,
<x_past>|month is the monthly average for a considered month.

The procedure is applied for each month to account for seasonality.

Methodology notes
-----------------
Originally the morphing approach originates from buildings physics, and is still in active
use for detailed simulation of heat transfer in buildings. See, for example:

J.B. Dias, G.C. da Graça, P.M.M. Soares (2020) Comparison of methodologies for
generation of future  building thermal energy simulation. Energy and Buildings,
vol. 206,2020, 109556, https://doi.org/10.1016/j.enbuild.2019.109556.
"""


import datetime as dt
import os

import atlite
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import configure_logging, create_logger
from shapely.geometry import LineString as Line
from shapely.geometry import Point

logger = create_logger(__name__)


def crop_cmip6(cmip6_xr, cutout_xr):
    """
    Crop the global-scale CMIP6 dataset to match the cutout extend.

    Parameters
    ----------
    cmip6_xr: xarray
        CMIP6 climate projections loaded as xarray
    cutout_xr : xarray of a cutout dataset
        Cutout dataset loaded as xarray
    """

    # spartial margin needed to avoid data losses during further interpolation
    d_pad = 1

    cmip6_region = cmip6_xr.sel(
        lat=slice(
            min(cutout_xr.coords["y"].values) - d_pad,
            max(cutout_xr.coords["y"].values) + d_pad,
        ),
        lon=slice(
            min(cutout_xr.coords["x"].values) - d_pad,
            max(cutout_xr.coords["x"].values) + d_pad,
        ),
    )
    return cmip6_region


def interpolate_cmip6_to_cutout_grid(cmip6_xr, cutout_xr):
    """
    Interpolate CMIP6 dataset to the cutout grid.

    Parameters
    ----------
    cmip6_xr: xarray
        CMIP6 climate projections loaded as xarray
    cutout_xr : xarray of a cutout dataset
        Cutout dataset loaded as xarray
    """

    # TODO read from the cutout file instead of hardcoding
    dx_new = 0.3

    newlon = np.arange(
        round(min(cutout_xr.coords["x"].values), 1),
        round(max(cutout_xr.coords["x"].values) + dx_new, 1),
        dx_new,
    )
    newlat = np.arange(
        round(min(cutout_xr.coords["y"].values), 1),
        round(max(cutout_xr.coords["y"].values) + dx_new, 1),
        dx_new,
    )
    cmip6_interp = cmip6_xr.interp(time=cmip6_xr.time, lat=newlat, lon=newlon)

    return cmip6_interp


# TODO fix years_window
def subset_by_time(cmip6_xr, month, year, years_window):
    """
    Filter CMIP6 dataset according to the requested year and the processed
    month.

    Parameters
    ----------
    cmip6_xr: xarray
        CMIP6 climate projections loaded as xarray
    month: integer
        A month value to be considered further for the morphing procedure
    year: integer
        The first year of the requested period in the future
    years_window: integer
        A width of the considered time period

    Example
    -------
    month=3, year=2050, years_window=30
    filters CMIP6 cmip6_xr dataset to calculate projection for March at 2050-2079
    """

    # TODO add warning in case there are no enough years
    cmip6_in_period = cmip6_xr.where(
        (
            (cmip6_xr["time.month"] == month)
            & (cmip6_xr["time.year"] >= year - years_window)
            & (cmip6_xr["time.year"] < year + years_window - 1)
        ),
        drop=True,
    )
    return cmip6_in_period


# TODO fix hardcoding in the parameter name
# TODO add functionality to customise CMIP6 ensemble
def calculate_proj_of_average(
    cmip6_xr, month, year0, year1, years_window, cmip6_param_name="t"
):
    """
    Calculate a change in a monthly average value for a climate parameter for
    the requested year periods and the processed month (a single value for each
    grid cell)

    Parameters
    ----------
    cmip6_xr: xarray
        CMIP6 climate projections loaded as xarray
    month: float
        A month value to be considered further for the morphing procedure
    year0: integer
        The first year of an earlier period in the future
    year1: integer
        The first year of a later period in the future
    years_window: integer
        Width of the considered time period

    Outputs
    -------
    dt_interp: xarray of the changes in the monthly averages for each grid cell
    Contains a single value for each grid cell.

    Example
    -------
    month=3, year0=2020, year0=2070, years_window=30
    calculates a multimodel change in the monthly average for
    March in 2020-2049 as compared with 2070-2099
    """

    cmip6_interp_year0 = subset_by_time(
        cmip6_xr, month, year=year0, years_window=years_window
    )
    cmip6_interp_year1 = subset_by_time(
        cmip6_xr, month=month, year=year1, years_window=years_window
    )
    dt_interp = cmip6_interp_year1[cmip6_param_name].mean("member").mean(
        "time"
    ) - cmip6_interp_year0[cmip6_param_name].mean("member").mean("time")
    return dt_interp


# TODO fix hardcoding in the parameter name
def build_projection_for_month(
    cutout_xr, dt_xr, month, cutout_param_name="temperature"
):
    """
    Calculate the projection for the cutout data applying the morphing approach
    using a dataset dt_xr for the changes in monthly aggregates for the
    processed month.

    Parameters
    ----------
    cutout_xr: xarray
        Cutout dataset loaded as xarray
    dt_xr: xarray
        Dataset of the changes in the monthly averages for each grid cell
    month: float
        A month value to be considered further for the morphing procedure

    Output
    -------
    cutout_xr: xarray of the projected changes for each grid cell for a processed month
    Contains a time series for each grid cell.
    """

    k_time = cutout_xr.time.dt.month.isin(month)

    for i in range(0, len(cutout_xr.y)):
        for j in range(0, len(cutout_xr.x)):
            cutout_xr[cutout_param_name][k_time, i, j] = np.add(
                cutout_xr[cutout_param_name][k_time, i, j],
                np.array([dt_xr[i, j].item()] * k_time.sum().item()),
            )

    return cutout_xr


def build_cutout_future(cutout_xr, cmip6_xr, months, year0, year1, years_window):
    """
    Calculate the projection for the cutout data applying the morphing approach
    for each month in the months.

    Parameters
    ----------
    cutout_xr: xarray
        Cutout dataset loaded as xarray
    cmip6_xr: xarray
        CMIP6 climate projections loaded as xarray
    months: list
        Months values to be used as a definition of a season to be considered
        for building projection time series
    year0: integer
        The first year of an earlier period in the future
    year1: integer
        The first year of a later period in the future
    years_window: integer
        Width of the considered time period

    Output
    -------
    cutout_xr: xarray
        Dataset on the projected changes for each grid cell for a processed month
        Contains a time series for each grid cell.
    """
    for k_month in [months]:
        dt_k = calculate_proj_of_average(
            cmip6_xr=cmip6_xr, month=k_month, year0=year0, year1=year1, years_window=5
        )

        build_projection_for_month(cutout_xr, dt_k, k_month)

    cutout_xr = cutout_xr.where(cutout_xr.time.dt.month.isin(months), drop=True)

    return cutout_xr


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "build_climate_projections",
            climate_scenario="ssp2-2.6",
            present_year=2020,
            future_year=2050,
            years_window=5,
            cutout="africa-2013-era5",
        )
    configure_logging(snakemake)

    climate_scenario = snakemake.params.climate_scenario
    present_year = snakemake.params.present_year
    future_year = snakemake.params.future_year
    years_window = snakemake.params.years_window
    cutout = snakemake.input.cutout
    cmip6 = snakemake.input.cmip6_avr

    snapshots = pd.date_range(freq="h", **snakemake.params.snapshots)
    season_in_focus = snapshots.month.unique().to_list()

    cmip6_xr = xr.open_dataset(cmip6)
    cutout_xr = xr.open_dataset(cutout)
    cmip6_region = crop_cmip6(cmip6_xr, cutout_xr)

    cmip6_region_interp = interpolate_cmip6_to_cutout_grid(
        cmip6_xr=cmip6_region, cutout_xr=cutout_xr
    )

    # -----------------------------------------------------------------
    #       to be replaced after debug
    # -----------------------------------------------------------------
    # graphical test of interpolation
    fig = plt.figure()
    cmip6_region_interp["t"].mean("time").mean("member").plot()
    fig.savefig("test_cmip6_interp.png", dpi=700)

    fig = plt.figure()
    (cutout_xr["temperature"] - 273.15).mean("time").plot(cmap="OrRd")
    fig.savefig("results/test_present.png", dpi=700)
    # -----------------------------------------------------------------

    # TODO Add a check for time match
    cutout_future = build_cutout_future(
        cutout_xr=cutout_xr,
        cmip6_xr=cmip6_region_interp,
        months=season_in_focus,
        year0=present_year,
        year1=future_year,
        years_window=years_window,
    )

    cutout_future.to_netcdf(snakemake.output[0])

    # -----------------------------------------------------------------
    #       to be replaced after debug
    # -----------------------------------------------------------------
    # graphical test of projection
    fig = plt.figure()
    (cutout_future["temperature"] - 273.15).mean("time").plot(cmap="OrRd")
    fig.savefig("results/test_cmip6_project.png", dpi=700)
