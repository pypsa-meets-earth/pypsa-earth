..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _structure:

##########################################
The data worflow
##########################################

**1. Grid topology data**
OpenStreetMap (OSM) data are being retrieved from the OSM server to be used to build the power grid topology model.
It's possible to check and validate OSM-extracted data directly:
- Snakemake scripts
- validation notebooks

**2. Climate data**
Currently we relay on ERA5 <link> data and `atlite` to process them.
Explain what is cutout and how it's used by generate_renewable_profiles

**3. General data**
There are a number GIS data which allow to build a realistic model. Currently we are using the following resources.

TODO Add credits, links and license informations

- *environmental*
GEBCO_2021_TID.nc
eez_v11.gpkg
copernicus/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif
landcover

- *economical*
costs.csv

- *technological*
custom_powerplants.csv
eia_hydro_annual_generation.csv

**4. Pre-calculated datasets**
There are some datasets which were prepared to ensure smooth run of the model. However, they may (and in some cases) must be replaced by custom ones. 

- geo-spatial data on location of protected and reserved areas

natura.tiff
Currently the pre-build file is calculated for Africa. 

- electricity demand profiles
data/ssp2-2.6/2030/era5_2013/Africa.nc



