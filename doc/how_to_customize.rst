..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _how_to_customize:

##########################################
How to customize?
##########################################

**1. Adjust the config**
The main parameters needed to customize the modeling process are defined in the configuration file `config.yaml`. Depending on the particular problem you are working on it could make sense to change a bit config settings. In particular, it's better to set the  as big as it's feasible when generating a cutout. A considered area can be narrowed anytime when building a specific model.

A proper match between the temporal and spatial parameters across the configuration file is essencial to build the model. Generally, if any misterious error message appears during the model set-up process there are chances that it can be resolved by a simple config check.

It could be helpful to keep in mind the following points:
1) the cutout name should be the same across the whole configuration file;
2) the countries of interest defined with `countries` list in the `config.yaml` should be covered by the cutout area;
3) the cutout time dimension, the weather year used for demand modeling and the actual snapshot should match.

**2. Load all the common data**
The `retrieve_databundle` rule makes the job.
Caution: mind the flags in the config

**3. Build the custom cutout**
What is cutout
Atlite approach is used 
Copernicis API is needed

Be careful when setting cutout coordinates and generally when `build_cutout` flag is on.

**4. Build a natura.tiff raster**
Is used to account for landuse restrictions on the protected and reserved nature areas
There is a pre-built one which is currently valid for Africa
It can be visualized
