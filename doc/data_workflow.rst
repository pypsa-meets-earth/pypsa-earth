..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _data_workflow:

##########################################
Data used by the model
##########################################

PyPSA-Earth has the flexible scope when a modeling domain can be selected from the whole Earth to any subregion. The modeling workflow includes data management workflow which is focused on open data collection, modification, prediction and validation. The developed data processing approaches aimed to provide model-ready data with temporal and spatial resolution as high as it's needed by a specific modeling task.

Here we'll look into architecture of the data workflow while practical hand-ons are given in the Tutorial section.

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
Data validation tips
===================================

The following validation points are worth keeping in mind when building your energy model:

1. Check the [power grid](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/network_validation.ipynb):
    - overall lines length;
    - general grid topology;
    - ensure that the general structure of the grid model is appropriate, playing with `tol` values and augmentation options if needed.
 
2. Compare the [installed capacity](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/capacity_validation.ipynb) values 

3. Validate the [power demand](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/demand_validation.ipynb) values and profile.

4. Check that [hydro](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/hydro_generation_validation.ipynb), [solar and wind](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/renewable_potential_validation.ipynb) potentials have reasonable values

5. Simulate the actual [energy mix](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/validation_nigeria.ipynb). Look for detailed explanations in https://arxiv.org/abs/2209.04663, section 5.1.

Data availability and quality usually is the biggest concern. Some useful hints on the real-world validation example can be found in the [Nigeria validation](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/validation_nigeria.ipynb) notebook.
