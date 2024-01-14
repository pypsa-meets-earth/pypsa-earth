.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _customization_basic1:

#######################
1. Basic customization
#######################

A good starting point to customize your model are settings of the default configuration file `config.default`. You may want to do a reserve copy of your current configuration file and then overwrite it by a default configuration:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % cp config.default.yaml config.yaml

The model can be adapted to include any country, multiple countries (e.g. Nigeria and Benin) or continents (currently `Africa` work as a whole continent) using `countries` argument:

.. code:: yaml

    countries: ["NG", "BJ"]

.. note::

    To build the model of regions outside of Africa, it is important to setup ``Copernicus Climate Data Store (CDS) API`` to build custom cutouts. 
    The same is true if the weather year other than 2013 is considered even for Africa. The detailed instructions for setting up the Copernicus API can be found in :ref:`customization_copernicus`.  

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