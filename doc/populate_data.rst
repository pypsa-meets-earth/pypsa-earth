.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

##########################################
Populate data
##########################################

Once the download and filtering process is complete, the clean data are then processed through 
a set of ``snakemake`` rules to identify the main data inputs for the network modelling.

The following list of rules apply:

- :mod:`build_cutout` prepares smaller weather data portions from `ERA5 <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`_ for cutout ``africa-2013-era5``
- :mod:`build_bus_regions` determines `Voronoi cells <https://en.wikipedia.org/wiki/Voronoi_diagram>`_ for all substations
- :mod:`build_powerplants` for today's thermal power plant capacities using `powerplantmatching <https://github.com/FRESNA/powerplantmatching>`_ allocating these to the closest substation for each powerplant
- :mod:`build_natura_raster` for rasterising `World Database on Protected Areas (WDPA) <https://www.protectedplanet.net/en/resources/wdpa-manual>`_
- :mod:`build_renewable_profiles` for the hourly capacity factors and installation potentials constrained by land-use in each substation's Voronoi cell for PV, onshore and offshore wind

Index:

.. toctree::
   :caption: Overview

   populate/build_cutout
   populate/build_bus_regions
   populate/build_powerplants
   populate/build_natura_raster
   populate/build_renewable_profiles
