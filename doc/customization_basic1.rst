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

To accurately model the temporal and spatial availability of renewables such as wind and solar energy, we process historical weather data using `atlite <https://atlite.readthedocs.io/en/latest/>`__ package.
In order to simulate a specific country, it is crucial to specify the `cutout region` within the atlite configuration in the ``config.yaml`` file:

.. code:: yaml

    atlite:
        nprocesses: 4
        cutouts:
            cutout-2013-era5-tutorial:
                module: era5
                dx: 0.3  # cutout resolution
                dy: 0.3  # cutout resolution
                # The cutout time is automatically set by the snapshot range.

Replace ``cutout-2013-era5-tutorial`` with the region of interest. For example, when simulating Kazakhstan, it should be updated to ``asia-2013-era5``.
Note please that a temporal dimension of the cutout should be consistent with the values set for `snapshots` parameter. A time range of the cutout is determined by the parameters set when building this cutout while the time resolution corresponds to those of the used climate archives. In case of ERA5 dataset used in PyPSA-Earth by default, hourly resolution is implied.
It is also possible to decide which weather data source should be used to calculate potentials and capacity factor time-series for each carrier.
For example, we may want to use the ERA-5 dataset for solar and not the default SARAH-2 dataset. Visit :ref:`config` page to get familiar with configuration details.

Finally, for countries beyond Africa, it is imperative to initially set ``build_cutout: true`` for the first run. This facilitates the construction of cutouts from weather data. Subsequently, it can be switched to false to avoid reconstructing the cutout. The same applies for other parameters in ``enable`` section.

.. code:: yaml

    enable:
        retrieve_databundle: true  #  Recommended 'true', for the first run. Otherwise data might be missing.
        retrieve_cost_data: true  # true: retrieves cost data from technology data and saves in resources/costs.csv, false: uses cost data in data/costs.csv
        download_osm_data: true  # If 'true', OpenStreetMap data will be downloaded for the above given countries
        build_natura_raster: true # If True, than an exclusion raster will be build
        build_cutout: false
        # If "build_cutout" : true, then environmental data is extracted according to `snapshots` date range and `countries`

To delve into the specifics of the provided configurations and explore additional settings, please refer to the :ref:`config` page.