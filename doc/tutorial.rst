.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _tutorial:


##########################################
Tutorial
##########################################

How to customise PyPSA-Earth?
=============================

A good starting point to customize your model are settings of the default configuration file `config.default`. You may want to do a reserve copy of your current configuration file and then overwrite it by a default configuration:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % cp config.default.yaml config.yaml

The model can be adapted to include any country, multiple countries (e.g. Nigeria and Benin) or continents (currently `Africa` work as a whole continent) using `countries` argument:

.. code:: yaml

    countries: ["NG", "BJ"]

Likewise, the example's temporal scope can be restricted (e.g. to 7 days):

.. code:: yaml

    snapshots:
        start: "2013-03-1"
        end: "2013-03-7"
        inclusive: "left" # end is not inclusive


Year-related parameters are also being used  when specifying `load_options`:

.. code:: yaml

    load_options:
      ssp: "ssp2-2.6"
      weather_year: 2013
      prediction_year: 2030
      scale: 1

The `weather_year` value corresponds to the weather data which was used to generate the electricity demand profiles for a selected area while `prediction_year` correspond to the point of a ssp trajectory. The available values for `weather_year` and `prediction_year` can be checked by looking into `pypsa-earth/data/ssp2-2.6` folder. Currently, there are pre-calculated demand data for 2011, 2013, 2018 weather years and for 2030, 2040, 2050, and 2100 scenario prediction years.

It is also possible to allow less or more carbon-dioxide emissions, while defining the current emissions.
It is possible to model a net-zero target by setting the `co2limit` to zero:

.. code:: yaml

    electricity:
        voltages: [220., 300., 380.]
        co2limit: 1.487e+9
        co2base: 1.487e+9

PyPSA-Earth can generate a database of existing conventional powerplants through open data sources.
It is possible to select which types of powerplants to be included:

.. code:: yaml

    extendable_carriers:
        Generator: [solar, onwind, offwind-ac, offwind-dc, OCGT]
        StorageUnit: [] # battery, H2
        Store: [battery, H2]
        Link: []  # H2 pipeline


To accurately model the temporal and spatial availability of renewables such as wind and solar energy, we rely on historical weather data.
It is advisable to adapt the required range of coordinates to the selection of countries.

.. code:: yaml

    atlite:
        nprocesses: 4
        cutouts:
                africa-2013-era5-tutorial:
                    module: era5
                    dx: 0.3  # cutout resolution
                    dy: 0.3  # cutout resolution
                    # The cutout time is automatically set by the snapshot range.


Note please that a temporal dimension of the cutout should be consistent with the values set for `snapshots` parameter. A time range of the cutout is determined by the parameters set when building this cutout while the time resolution corresponds to those of the used climate archives. In case of ERA5 dataset used in PyPSA-Earth by default, hourly resolution is implied.

It is also possible to decide which weather data source should be used to calculate potentials and capacity factor time-series for each carrier.
For example, we may want to use the ERA-5 dataset for solar and not the default SARAH-2 dataset.

.. code:: yaml

    africa-2013-era5-tutorial:
        module: era5

.. code:: yaml

    solar:
        cutout: africa-2013-era5-tutorial

Finally, it is possible to pick a solver. For instance, this tutorial uses the open-source solver glpk and does not rely
on the commercial solvers such as Gurobi or CPLEX (for which free academic licenses are available).

.. code:: yaml

    solver:
        name: glpk


Be mindful that we only noted major changes to the provided default configuration that is comprehensibly documented in :ref:`config`.
There are many more configuration options beyond what is adapted for the tutorial!

How to execute different parts of the workflow?
===============================================

Snakemake is a workflow management tool inherited by PyPSA-Earth from PyPSA-Eur.
Snakemake decomposes a large software process into a set of subtasks, or ’rules’, that are automatically chained to obtain the desired output.

.. note::

  ``Snakemake``, which is one of the major dependencies, will be automatically installed in the environment pypsa-earth, thereby there is no need to install it manually.

The snakemake included in the conda environment pypsa-earth can be used to execute any custom rule with the following command:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake < your custom rule >  

Starting with essential usability features, the implemented PyPSA-Earth `Snakemake procedure <https://github.com/pypsa-meets-earth/pypsa-earth/blob/main/Snakefile>`_ that allows to flexibly execute the entire workflow with various options without writing a single line of python code. For instance, you can model the world energy system or any subset of countries only using the required data. Wildcards, which are special generic keys that can assume multiple values depending on the configuration options, help to execute large workflows with parameter sweeps and various options.

You can execute some parts of the workflow in case you are interested in some specific it's parts.
E.g. power grid topology may be extracted and cleaned with the following command which refers to the script name: 

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 clean_osm_data

Solar profile for the requested area may be calculated using the output name:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 resources/renewable_profiles/profile_solar.nc 


How to use PyPSA-Earth for your energy problem?
===============================================

PyPSA-Earth mostly relies on the :ref:`global datasets <data_workflow>` and can be tailored to represent any part of the world in a few steps. The following procedure is recommended.

1. Adjust the model configuration
---------------------------------

The main parameters needed to customize the inputs for your national-specific data are defined in the :ref:`configuration <config>` file `config.yaml`. The configuration settings should be adjusted according to a particular problem you are intended to model. The main regional-dependent parameters are:

* `countries` parameter which defines a set of the countries to be included into the model;

* `cutouts` and `cutout` parameters which refer to a name of the climate data archive (so called `cutout <https://atlite.readthedocs.io/en/latest/ref_api.html#cutout>`_) to be used for calculation of the renewable potential.

Apart of that, it's worth to check that there is a proper match between the temporal and spatial parameters across the configuration file as it is essential to build the model properly. Generally, if there are any mysterious error message appearing during the first model run, there are chances that it can be resolved by a simple config check.

It could be helpful to keep in mind the following points:

1. the cutout name should be the same across the whole configuration file (there are several entries, one under under `atlite` and some under each of the `renewable` parameters);

2. the countries of interest defined with `countries` list in the `config.yaml` should be covered by the cutout area;

3. the cutout time dimension, the weather year used for demand modeling and the actual snapshot should match.

2. Build the custom cutout
--------------------------

The cutout is the main concept of climate data management in PyPSA ecosystem introduced in `atlite <https://atlite.readthedocs.io/en/latest/>`_ package. The cutout is an archive containing a spatio-temporal subset of one or more topology and weather datasets. Since such datasets are typically global and span multiple decades, the Cutout class allows atlite to reduce the scope to a more manageable size. More details about the climate data processing concepts are contained in `JOSS paper <https://joss.theoj.org/papers/10.21105/joss.03294>`_.

Generally, the spatial and time resolution of the cutout data is determined by parameters of an underlying dataset. That is 30 km x 30 km grid and houtly resolution for ERA5 recommended for usage in PyPSA-Earth.

The pre-built cutout for Africa is available for 2013 year and can be loaded directly from zenodo through the rule `retrieve_cutout`. There is also a smaller cutout for Africa built for a two-weeks time span; it is automatically downloaded when retrieving common data with `retrieve_databundle_light`.

.. note::
    Skip this recommendation if the region of your interest is within Africa and you are fine with the 2013 weather year

In case you are interested in other parts of the world you have to generate a cutout yourself using the `build_cutouts` rule. To run it you will need to: 

1. be registered on  the `Copernicus Climate Data Store <https://cds.climate.copernicus.eu>`_;

2. install `cdsapi` package  (can be installed with `pip`);

3. setup your CDS API key as described `on their website <https://cds.climate.copernicus.eu/api-how-to>`_.

These steps are required to use CDS API which allows an automatic file download while executing `build_cutouts` rule.

The `build_cutout` flag should be set `true` to generate the cutout. After the cutout is ready, it's recommended to set `build_cutout` to `false` to avoid overwriting the existing cutout by accident. The `snapshots` values set when generating the cutout, will determine the temporal parameters of the cutout. Accessible years which can be used to build a cutout depend on ERA5 data availability. `ERA5 page <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`_ explains that the data is available from 1950 and updated continuously with about 3 month delay while the data on 1950-1978 should be treated as preliminary as that is a rather recent development.

After the first run, if you don't change country and don't need to increase a considered time span wider than the one you created the cutout with, you may set to false both `retrieve_databundle` and `build_cutout`.

Spatial extent
^^^^^^^^^^^^^^

Normally cutout extent is calculated from the shape of the requested region defined by the `countries` parameter in the configuration file `config.yaml`. It could make sense to set the countries list as big as it's feasible when generating a cutout. A considered area can be narrowed anytime when building a specific model by adjusting content of the `countries` list.

There is also option to set the cutout extent specifying `x` and `y` values directly. However, these values will overwrite values extracted from the countries shape. Which means that nothing prevents `build_cutout` to extract data which has no relation to the requested countries. Please use direct definition of `x` and `y` only if you really understand what and why you are doing.

Temporal extent
^^^^^^^^^^^^^^^

If you create the cutout for a certain year (let's say 2013) and want to run scenarios for a subset of this year, you don't need to rerun the `build_cutout` as the cutout still contains all the hours of 2013. The workflow will automatically subset the cutout archive to extract data for the particular timeframe of interest. If you instead you want to run the 2014 scenario, then rerun `build_cutout` is needed. 

In case you need model a number of years, a convenient approach may be to create the cutout for the whole period under interest (e.g. 2013-2015) so that you don't need to build any additional cutouts. Note, however, that the disk requirements increase in this case.

3. Check natura.tiff raster
-----------------------------

A raster file `natura.tiff` is used to store shapes of the protected and reserved nature areas. Such landuse restrictions can be taking into account when calculating the renewable potential with `build_renewable_profiles` by switching-on `natura` option:

.. code:: bash

    natura: true

A pre-built `natura.tiff` is loaded along with other data needed to run a model with `retrieve_databundle_light` rule. This raster file contains data on the on protected areas around the world where areas no (renewable) assets can be installed. The `natura.tiff` raster has now global coverage so you don't need to create it locally.

How to validate?
================

.. TODO add a list of actions needed to do the validation

To validate the data obtained with PyPSA-Earth, we recommend to go through the procedure here detailed. An exampled of the validation procedure is available in the `Nigeria validation <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/validation_nigeria.ipynb>`_ notebook. Public information on the power system of Nigeria are compared to those obtained from the PyPSA-Earth model.

Simulation procedure
--------------------

It may be recommended to check the following quantities the validation:

#. inputs used by the model:

    #. network characteristics;

    #. substations;

    #. installed generation by type;

#. outputs of the simulation:

    #. demand;

    #. energy mix.

Where to look for reference data
--------------------------------
 
Data availability for many parts of the world is still quite limited. Usually the best sources to compare with are regional data hubs. There is also a collection of harmonized datasets curated by the international organisations. A non-exhaustive list of helpful sources:

* `World Bank <https://energydata.info/>`_;

* International Renewable Energy Agency `IRENA <https://pxweb.irena.org/pxweb/en/IRENASTAT/IRENASTAT__Power%20Capacity%20and%20Generation/ELECCAP_2022_cycle2.px/>`_;

* International Energy Agency `IEA <https://www.iea.org/data-and-statistics>`_;

* `BP <https://www.bp.com/en/global/corporate/energy-economics/statistical-review-of-world-energy.html>`_ Statistical Review of World Energy;

* International Energy Agency `IEA <https://www.iea.org/data-and-statistics>`_;

* `Ember <https://ember-climate.org/data/data-explorer/>`_ Data Explorer.


Advanced validation examples
----------------------------

The following validation notebooks are worth a look when validating your energy model:

1. A detailed `network validation <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/network_validation.ipynb>`_.
 
2. Analysis of `the installed capacity <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/capacity_validation.ipynb>`_ for the considered area. 

3. Validation of `the power demand <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/demand_validation.ipynb>`_ values and profile.

4. Validation of `hydro <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/hydro_generation_validation.ipynb>`_, `solar and wind <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/renewable_potential_validation.ipynb>`_ potentials.


.. include:: ./how_to_docs.rst
