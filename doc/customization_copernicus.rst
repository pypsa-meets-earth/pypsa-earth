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