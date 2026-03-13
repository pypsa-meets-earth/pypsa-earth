# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import os

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr

# cutout to be checked
cutout_fl_path = "cutouts/cutout-2013-era5.nc"

data = xr.open_dataset(os.path.join(cutout_fl_path))

## Basic plotting

### Quickly cross-check functionality

data["wnd100m"].mean("time").plot.contourf(cmap=plt.cm.Blues)
data["temperature"].mean("time").plot.contourf(cmap=plt.cm.autumn)


## Advanced plotting

### Make sure the spatial distribution makes sense for a region of interest

import cartopy

### Plot wind
fig, ax = plt.subplots(figsize=(8, 4), subplot_kw={"projection": ccrs.PlateCarree()})
data["wnd100m"].mean("time").plot.contourf(cmap=plt.cm.Blues)
ax.add_feature(cartopy.feature.BORDERS, color="gray", linewidth=0.5)
ax.coastlines()
plt.savefig("cutout_wind.jpeg")

### Plot solar
fig, ax = plt.subplots(figsize=(8, 4), subplot_kw={"projection": ccrs.PlateCarree()})
data["influx_direct"].mean("time").plot.contourf(cmap=plt.cm.Oranges)
ax.add_feature(cartopy.feature.BORDERS, color="gray", linewidth=0.5)
ax.coastlines()

### Plot temperature
fig, ax = plt.subplots(figsize=(8, 4), subplot_kw={"projection": ccrs.PlateCarree()})
data["temperature"].mean("time").plot.contourf(cmap=plt.cm.summer)
ax.add_feature(cartopy.feature.BORDERS, color="gray", linewidth=0.5)
ax.coastlines()
plt.savefig("cutout_tas.jpeg")
