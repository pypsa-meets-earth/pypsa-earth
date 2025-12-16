# Download & Filter

The download and filtering process of the PyPSA-Earth energy system model consists of a group of `snakemake` rules which are briefly outlined and explained in detail in the sections below.

Not all data dependencies are shipped with the git repository. Instead we provide separate data bundles which are loaded automatically when running `solve_all_networks` rule when `retrieve_databundle_light` flag in the configuration file is on.

## Rules

- **[retrieve_databundle_light](retrieve-databundle-light.md)** - Enables a simplified approach to download the main databundle of raw data.

- **[build_shapes](build-shapes.md)** - Automatically downloads administrative country shapes from the [GADM dataset](https://gadm.org/) and generates GeoJSON files with shapes of the countries, exclusive economic zones and administrative zones at the desired resolution (e.g. region, district, etc.).
