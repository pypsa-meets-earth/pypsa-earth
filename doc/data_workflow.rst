..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _data_workflow:

##########################################
Data used by the model
##########################################

PyPSA-Earth has the flexible scope when a modeling domain can be selected from the whole Earth to any subregion. The modeling workflow includes data management workflow which is focused on open data collection, modification, prediction and validation. The developed data processing approaches aimed to provide model-ready data with temporal and spatial resolution as high as it's needed by a specific modeling task.

.. _data_management_strategy:

Data management strategy
===================================

Here we'll look into architecture of the data workflow while practical hand-ons are given in the Tutorial section.

**1. Grid topology data**
OpenStreetMap `OSM <https://www.openstreetmap.org/>`_ data are used to build power grid topology model. OSM is the biggest crowd-sourced collection of geographic information, which is daily updated and includes geolocation references. The OSM data are being loaded by `download_osm_data` and cleaned by `clean_osm_data` rules, respectively.

**2. Climate data**
The climate data processing is provided by `atlite <https://atlite.readthedocs.io/en/latest/>`_ package. It extracts all the required whether and climate data to generate the time series of renewable potential by g`enerate_renewable_profiles` rule.

The main data source on climate variables is `ERA5 reanalysis <https://rmets.onlinelibrary.wiley.com/doi/10.1002/qj.3803>`_.

**3. General data**
There are a number datasets applied in PyPSA-Earth to build a realistic model. Original datasets are stored in the `./pypsa-earth/data/` folder.

Currently we are using the following resources.

- *environmental*

**copernicus** contains the raw data on the land covering as available from the Copernicus database. It is used in the build_renewable_profiles rule to quantify what are the land regions available for the installation of renewable resources, e.g. renewable assets may not be installed on arable land

**eez** is the dataset of the Exclusive Economic Zones (EEZ) available from Marine Regions. This file is used in the rule build_shapes to identify the marine region by country and provide shapes of the maritime regions to be possibly used to estimate off-shore renewable potential, for example.

**gebco** gridded bathymetric data which can be translated into depths and shapes of underwater terrain. These data are used in the `build_renewable_profiles` rule. [GEBCO](https://www.gebco.net/) stands for General Bathymetric Chart of the Oceans. It's curated by a non-profit making organisation which relies largely on the voluntary contributions of an enthusiastic international team of geoscientists and hydrographers.

**hydrobasins** datasets on watershed boundaries and basins, as available from HydroBASINS. These data are used to estimate the hydropower generation in the `build_renewable_profiles` rule.

**landcover** describes the shapes of world protected areas that are needed to identify in what areas no (renewable) assets can be installed. Currently are used to generate a `natura.tiff` raster. Will be deprecated once the global `natura.tiff` will be available.

**osm** are raw [OpenStreetMap](https://www.openstreetmap.org/) data. Are being loaded behind the scene when running the `download_osm_data` rule.

The raw OSM data in [.pbf](https://wiki.openstreetmap.org/wiki/PBF_Format) format are stored in the folder `./pypsa-earth/data/osm/{region}/pbf/`. Here `region` denotes a continent, e.g. Africa, or a macro region, e.g. Central America, where the countries of interest belong. The pbf-files contain the entire OSM data for the country; the specific network information related to generators, substations, lines and cables are extracted, cleaned and writen as a geojsons in the folder `./pypsa-earth/data/osm/{region}/Elements/`. All network data (generators, substations, lines and cables) for each country are stored as geojson files.

The cleaned OSM network data that are the output of the `osm_data_cleaning` rule, which process the raw OSM data to obtain cleaned datasets of all the network assets, namely generators, substations, lines and cables. These data are stored in `/resources/osm/` folder.

- *economical*

**costs.csv**
csv file containing the defaulf costs of the technologies along with their typical lifetime and efficiency values. The dataset is intended to give a starting point for running the model while regional adajustments may be needed. 

**gadm** it contains data of the shapes of administrative zones by country (e.g. regions, districts, provinces, ...), depending on the level of resolution desired by the configuration file. The data in this folder are automatically populated by the build_shapes rule that download such data from the gadm website

**GDP** raster dataset of the Gross Domestic Product (GDP) by arcs of the world, as available from DRYAD

**WorldPop** raster dataset of the population by arc as automatically by build_shapes rule from WorldPop

- *technological*
**custom_powerplants.csv**
**eia_hydro_annual_generation.csv**
**hydro_capacities.csv**

**4. Pre-calculated datasets**
There are some datasets which were prepared to ensure smooth run of the model. However, they may (and in some cases) must be replaced by custom ones. 

- geo-spatial data on location of protected and reserved areas

natura.tiff
Currently the pre-build file is calculated for Africa. 

- electricity demand profiles
data/ssp2-2.6/2030/era5_2013/Africa.nc

.. .. _data_validation_tips:

.. Data validation tips
.. ===================================

.. The following validation points are worth keeping in mind when building your energy model:

.. 1. Check the [power grid](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/network_validation.ipynb):
..     - overall lines length;
..     - general grid topology;
..     - ensure that the general structure of the grid model is appropriate, playing with `tol` values and augmentation options if needed.
 
.. 2. Compare the [installed capacity](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/capacity_validation.ipynb) values 

.. 3. Validate the [power demand](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/demand_validation.ipynb) values and profile.

.. 4. Check that [hydro](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/hydro_generation_validation.ipynb), [solar and wind](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/renewable_potential_validation.ipynb) potentials have reasonable values

.. 5. Simulate the actual [energy mix](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/validation_nigeria.ipynb). Look for detailed explanations in https://arxiv.org/abs/2209.04663, section 5.1.

.. Data availability and quality usually is the biggest concern. Some useful hints on the real-world validation example can be found in the [Nigeria validation](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/validation_nigeria.ipynb) notebook.
