..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _how_to_customize:

##########################################
How to customize?
##########################################

PyPSA-Earth can be tailored to represent any part of the world quite easily. The following procedure may be recommended.

**1. Adjust the config**
The main parameters needed to customize the modeling process are defined in the configuration file `config.yaml`. The configuration settings should be adjusted according to a particular problem you are intended to model. The main regional-dependent parameters are:
- `countries` parameter which defines a set of the countries to be included into the model;
- `cutouts` and `cutout` parameters which refer to a name of the climate data archive (so called *cutout*) to be used for calculation of the renewable potential.

Apart of that, it's worth to check that there is a proper match between the temporal and spatial parameters across the configuration file as it is essencial to build the model properly. Generally, if there are any misterious error message appearing during the first model run, there are chances that it can be resolved by a simple config check.

It could be helpful to keep in mind the following points:
1) the cutout name should be the same across the whole configuration file (there are several entries, one under under `atlite` and some under each of the `renewable` parameters);
2) the countries of interest defined with `countries` list in the `config.yaml` should be covered by the cutout area;
3) the cutout time dimension, the weather year used for demand modeling and the actual snapshot should match.

**2. Load all the common data**

PyPSA-Earth relies on a number of datasets introduced in the (:ref:`_data_workflow`) and (:ref:`_introduction#License`). Automated data retrival is implemented by the retrieve* rules (:ref:`rules_overview/download_and_filter/retrieve`).

Mostly these data are global and should be loaded only once when installing the model. The exception is the cutout which requires of about ~200 Gb in it's global version. The cutout which is supplied with the model could be loaded by running `retrieve_databundle_light` rule, corresponds to Africa. If you are interested in other part of the world a custom cutout should be built. 

**3. Build the custom cutout**
What is cutout
Atlite approach is used 
Copernicis API is needed

Be careful when setting cutout coordinates and generally when `build_cutout` flag is on.

**4. Build a natura.tiff raster**
Is used to account for landuse restrictions on the protected and reserved nature areas
There is a pre-built one which is currently valid for Africa
It can be visualized
