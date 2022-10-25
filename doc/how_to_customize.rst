..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _how_to_customize:

##########################################
How to customize?
##########################################

PyPSA-Earth can be tailored to represent any part of the world quite easily. The following procedure may be recommended.

**1. Adjust the config**
The main parameters needed to customize the inputs for your national-specific data are defined in the configuration file `config.yaml`. The configuration settings should be adjusted according to a particular problem you are intended to model. The main regional-dependent parameters are:
- `countries` parameter which defines a set of the countries to be included into the model;
- `cutouts` and `cutout` parameters which refer to a name of the climate data archive (so called *cutout*) to be used for calculation of the renewable potential.

Apart of that, it's worth to check that there is a proper match between the temporal and spatial parameters across the configuration file as it is essential to build the model properly. Generally, if there are any mysterious error message appearing during the first model run, there are chances that it can be resolved by a simple config check.

It could be helpful to keep in mind the following points:
1) the cutout name should be the same across the whole configuration file (there are several entries, one under under `atlite` and some under each of the `renewable` parameters);
2) the countries of interest defined with `countries` list in the `config.yaml` should be covered by the cutout area;
3) the cutout time dimension, the weather year used for demand modeling and the actual snapshot should match.

**2. Load all the common data**

PyPSA-Earth relies on a number of datasets introduced in the :ref:`_data_workflow` and :ref:`_introduction#License`. Automated data retrieval is implemented by the retrieve* rules :ref:`rules_overview/download_and_filter/retrieve`.

Mostly these data are global and should be loaded only once when installing the model. The exception is the cutout which requires of about ~200 Gb in it's global version. The cutout which is supplied with the model could be loaded by running `retrieve_databundle_light` rule, corresponds to Africa. If you are interested in other part of the world a custom cutout should be built. However, in the future we plan to provide a general default global climate inputs.

**3. Build the custom cutout**
The cutout is the main concept of climate data management in PyPSA ecosystem introduced in `atlite <https://atlite.readthedocs.io/en/latest/>`_ package. The cutout is an archive containing a spatio-temporal subset of one or more topology and weather datasets. Since such datasets are typically global and span multiple decades, the Cutout class allows atlite to reduce the scope to a more manageable size. More details about the climate data processing concepts are contained in `JOSS paper <https://joss.theoj.org/papers/10.21105/joss.03294>`_.

The pre-built cutout for Africa is available for 2013 year and can be loaded directly from zenodo through the rule `retrieve_cutout`. There is also a smaller cutout for Africa built for a two-weeks time span; it is automatically downloaded when retrieving common data with `retrieve_databundle_light`.

In case you are interested in other parts of the world you can generate a cutout yourself using the `build_cutouts` rule. To run it you will need to 
1) be registered on  the `Copernicus Climate Data Store <https://cds.climate.copernicus.eu>`_; 
2) install `cdsapi` package  (can be installed with `pip`);
3) setup your CDS API key as described `on their website <https://cds.climate.copernicus.eu/api-how-to>`_.

These steps are required to use CDS API which allows an automatic file download while executing `build_cutouts` rule.

Normally cutout extent is calculated from the shape of the requested region defined by the `countries` parameter in the configuration file `config.yaml`. It could make sense to set the countries list as big as it's feasible when generating a cutout. A considered area can be narrowed anytime when building a specific model by adjusting content of the `countries` list.

There is also option to set the cutout extent specifying `x` and `y` values directly. However, these values will overwrite values extracted from the countries shape. Which means that nothing prevents `build_cutout` to extract data which has no relation to the requested countries. Please use direct definition of `x` and `y` only if you really understand what and why you are doing.

The `build_cutout` flag should be set `true` to generate the cutout. After the cutout is ready, it's recommended to set `build_cutout` to `false` to avoid overwriting the existing cutout by accident.

**4. Build a natura.tiff raster**
A raster file `natura.tiff` is used to store shapes of the protected and reserved nature areas. Such landuse restrictions can be taking into account when calculating the renewable potential with `build_renewable_profiles`.

A pre-built `natura.tiff` is loaded along with other data needed to run a model with `retrieve_databundle_light` rule. Currently this raster is valid for Africa, global `natura.tiff` raster is under development. You may generate the `natura.tiff` for a region of interest using `build_natura_raster` rule which aggregates data on protected areas along the cutout extent.
