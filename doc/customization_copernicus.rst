.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _customization_copernicus:

####################
Setup Copernicus API
####################


To build the model of regions outside of Africa, it is important to access to `Copernicus Climate Data Store <https://cds.climate.copernicus.eu>`_ to build custom cutouts.
The same is true if the weather year other than 2013 is considered for the region of interest in Africa.

.. note::

    Skip this recommendation if the region of your interest is within Africa and you are fine with the 2013 weather year

Steps to get access to Copernicus database:

1. Register to  the `Copernicus Climate Data Store <https://cds.climate.copernicus.eu>`_;
2. Install `cdsapi` package  (can be installed with `pip`);
3. Setup your CDS API key as described `on their website <https://cds.climate.copernicus.eu/api-how-to>`_.

These steps are required to use CDS API which allows an automatic file download while executing `build_cutouts` rule.

The `build_cutout` flag should be set `true` to generate the cutout. After the cutout is ready, it's recommended to set `build_cutout` to `false` to avoid overwriting the existing cutout by accident. The `snapshots` values set when generating the cutout, will determine the temporal parameters of the cutout. Accessible years which can be used to build a cutout depend on ERA5 data availability. `ERA5 page <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`_ explains that the data is available from 1950 and updated continuously with about 3 month delay while the data on 1950-1978 should be treated as preliminary as that is a rather recent development.

After the first run, if you don't change country and don't need to increase a considered time span wider than the one you created the cutout with, you may set both `retrieve_databundle` and `build_cutout` to false.
