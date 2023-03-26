.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

##########################################
Download and filter data
##########################################

The download and filtering process of the PyPSA-Earth energy system model consists of a group of ``snakemake`` rules which are briefly outlined and explained in detail in the sections below.

Not all data dependencies are shipped with the git repository. Instead we provide separate data bundles which is loaded automatically when running `solve_all_networks` rule when `retrieve_databundle` flag in the configuration file is on. 

- :mod:`retrieve` enables a simplified approach to download the main databundle of raw data
- :mod:`build_shapes` automatically downloads administrative country shapes from the
  `GADM dataset <https://gadm.org/>`_ and generates GeoJSON files with shapes of the countries, exclusive economic zones and administrative zones at the desired resolution (e.g. region, district, etc.).


.. toctree::
   :caption: Overview

   download_and_filter/retrieve
   download_and_filter/build_shapes
