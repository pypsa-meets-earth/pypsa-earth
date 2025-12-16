# Populate Data

Once the download and filtering process is complete, the clean data are then processed through a set of `snakemake` rules to identify the main data inputs for the network modelling.

## Rules

- **[build_cutout](build-cutout.md)** - Prepares smaller weather data portions from [ERA5](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5) for renewable resource calculations.

- **[build_bus_regions](build-bus-regions.md)** - Determines [Voronoi cells](https://en.wikipedia.org/wiki/Voronoi_diagram) for all substations.

- **[build_powerplants](build-powerplants.md)** - Builds today's thermal power plant capacities using [powerplantmatching](https://github.com/FRESNA/powerplantmatching), allocating these to the closest substation for each powerplant.

- **[build_natura_raster](build-natura-raster.md)** - Rasterises [World Database on Protected Areas (WDPA)](https://www.protectedplanet.net/en/resources/wdpa-manual) data.

- **[build_renewable_profiles](build-renewable-profiles.md)** - Builds the hourly capacity factors and installation potentials constrained by land-use in each substation's Voronoi cell for PV, onshore and offshore wind.

- **[build_demand_profiles](build-demand-profiles.md)** - Builds the hourly demand profiles for each substation.
