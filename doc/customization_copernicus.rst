.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _customization_copernicus:

#######################
5. Setup Copernicus API
#######################


To build custom cutouts, it is important to access to `Copernicus Climate Data Store <https://cds.climate.copernicus.eu>`__.

.. note::

    Skip this recommendation if you are fine with the 2013 weather year.

Steps to get access to Copernicus database:

1. Register to  the `Copernicus Climate Data Store <https://cds.climate.copernicus.eu>`_;
2. Setup your CDS API key as described `on their website <https://cds.climate.copernicus.eu/how-to-api>`_.

These steps are required to use CDS API which allows an automatic file download while executing `build_cutouts` rule.

The `build_cutout` flag should be set `true` to generate the cutout. After the cutout is ready, it's recommended to set `build_cutout` to `false` to avoid overwriting the existing cutout by accident. The `snapshots` values set when generating the cutout, will determine the temporal parameters of the cutout. Accessible years which can be used to build a cutout depend on ERA5 data availability. `ERA5 page <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`_ explains that the data is available from 1950 and updated continuously with about 3 month delay while the data on 1950-1978 should be treated as preliminary as that is a rather recent development.

.. note::

    Building a cutout may require a significant amount of time and storage space. Continental cutouts, such as those for Asia, South America, and Africa, typically require around 20 GB of storage space, while cutouts for individual countries or small regions may occupy approximately 1-5 GB. The process of building a cutout can take between 2 to 3 hours.

After the first run, if you don't change country and don't need to increase a considered time span wider than the one you created the cutout with, you may set both `retrieve_databundle` and `build_cutout` to false.
