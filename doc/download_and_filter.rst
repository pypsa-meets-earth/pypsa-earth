..
  SPDX-FileCopyrightText: 2019-2020 The PyPSA-Africa Authors,
  adapted from PyPSA-Eur

  SPDX-License-Identifier: CC-BY-4.0

##########################################
Download and filter data
##########################################

The download and filtering process of the PyPSA-Africa energy system model consists of a 
group of ``snakemake`` rules which are briefly outlined and explained in detail in the
sections below.

Not all data dependencies are shipped with the git repository.
Instead we provide separate data bundles which can be obtained
using the ``retrieve*`` rules (:ref:`data`).
Having downloaded the necessary data,

- :mod:`retrieve` enables a simplified approach to download the main databundle of raw data
  needed for the execution
- :mod:`build_shapes` automatically downloads administrative country shapes from the
  `GADM dataset <https://gadm.org/>`_ and generates GeoJSON files with shapes of the countries,
  exclusive economic zones and administrative zones at the desired resolution
  (e.g. region, district, etc.).

Index:

.. toctree::
   :caption: Overview

   download_and_filter/retrieve
   download_and_filter/build_shapes
