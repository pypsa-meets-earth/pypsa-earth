


Configuration


PyPSA-Earth imports the configuration options originally developed in [PyPSA-Eur](https://pypsa-eur.readthedocs.io/en/latest/index.html) and here reported and adapted.
The options here described are collected in a [`config.yaml` file located in the root directory.
Users should copy the provided default configuration (`config.default.yaml`) and amend
their own modifications and assumptions in the user-specific configuration file (`config.yaml``);
confer installation instructions at [installation].


  Credits to PyPSA-Eur developers for the initial drafting of the configuration documentation here reported



# Top-level configuration


   :language: yaml
   :start-at: version:
   :end-at:  build_cutout:



   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/toplevel.csv



# `run`

It is common conduct to analyse energy system optimisation models for **multiple scenarios** for a variety of reasons,
e.g. assessing their sensitivity towards changing the temporal and/or geographical resolution or investigating how
investment changes as more ambitious greenhouse-gas emission reduction targets are applied.

The `run` section is used for running and storing scenarios with different configurations which are not covered by [wildcards]. It determines the path at which resources, networks and results are stored. Therefore the user can run different configurations within the same directory. If a run with a non-empty name should use cutouts shared across runs, set `shared_cutouts` to `true`.


   :language: yaml
   :start-at: run:
   :end-at: shared_cutouts:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/run.csv



# `scenario`

The `scenario` section is an extraordinary section of the config file
that is strongly connected to the [wildcards[ and is designed to
facilitate running multiple scenarios through a single command

``bash
snakemake -j 1 solve_all_networks

``

For each wildcard, a **list of values** is provided. The rule `solve_all_networks` will trigger the rules for creating `results/networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc` for **all combinations** of the provided wildcard values as defined by Python's `itertools.product(...)](https://docs.python.org/2/library/itertools.html#itertools.product) function that snakemake's [expand(...) function](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#targets) uses.

An exemplary dependency graph (starting from the simplification rules) then looks like this:

![Image](assets/images/scenarios.png)


   :language: yaml
   :start-at: scenario:
   :end-at: opts:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/scenario.csv



# [`snapshots``

Specifies the temporal range for the historical weather data, which is used to build the energy system model. It uses arguments to `pandas.date_range](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.date_range.html). The date range must be in the past (before 2022). A well-tested year is 2013.


   :language: yaml
   :start-at: snapshots:
   :end-at: inclusive:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/snapshots.csv



# [`crs`

Defines the coordinate reference systems (crs).


   :language: yaml
   :start-at: crs:
   :end-at: area_crs:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/crs.csv


# `augmented_line_connection``

If enabled, it increases the connectivity of the network. It makes the network graph `k-edge-connected](https://en.wikipedia.org/wiki/K-edge-connected_graph), i.e.,
if fewer than k edges are removed, the network graph stays connected. It uses the [k-edge-augmentation](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.connectivity.edge_augmentation.k_edge_augmentation.html#networkx.algorithms.connectivity.edge_augmentation.k_edge_augmentation)
algorithm from the [NetworkX](https://networkx.org/documentation/stable/index.html) Python package.


   :language: yaml
   :start-at: augmented_line_connection:
   :end-at: min_DC_length:


   :header-rows: 1
   :widths: 25,10,22,27
   :file: configtables/augmented_line_connection.csv



# [`cluster_options`

Specifies the options to simplify and cluster the network. This is done in two stages, first using the rule `simplify_network` and then using the rule `cluster_network``. For more details on this process, see the `PyPSA-Earth paper](https://www.sciencedirect.com/science/article/pii/S0306261923004609), section 3.7.


   :language: yaml
   :start-at: cluster_options:
   :end-at: efficiency:


   :header-rows: 1
   :widths: 25,10,22,27
   :file: configtables/cluster_options.csv



# [`build_shape_options`

Specifies the options to build the shapes in which the region of interest (`countries`) is divided.


   :language: yaml
   :start-at: build_shape_options:
   :end-at: contended_flag:


   :header-rows: 1
   :widths: 25,10,22,27
   :file: configtables/build_shape_options.csv



# `subregion`

If enabled, this option allows a region of interest (`countries`) to be redefined into subregions,
which can be activated at various stages of the workflow. Currently, it is used in `simplify_network` and `cluster_network` rule.


   :language: yaml
   :start-at: subregion:
   :end-at: path_custom_shapes:


   :header-rows: 1
   :widths: 25,10,22,27
   :file: configtables/subregion.csv

The names of subregions are arbitrary. Its sizes are determined by how many GADM IDs that are included in the list.
A single country can be divided into multiple subregions, and a single subregion can include GADM IDs from multiple countries.
If the same GADM ID appears in different subregions, the first subregion listed will take precedence over that region.
The remaining GADM IDs that are not listed will be merged back to form the remaining parts of their respective countries.
For example, consider the Central District of Botswana, which has a GADM ID of `BW.3`. To separate this district from the rest of the country, you can select:


   :language: yaml
   :start-at: subregion:
   :end-at: Central: 0.3

There are several formats for GADM IDs depending on the version, so before using this feature, please review the `resources/shapes/gadm_shape.geojson` file which can be created using the command:

``bash
snakemake -j 1 build_shapes

``


   The rule `build_shapes`` currently use `Version 4.1 ](https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/) for their GADM data. This may change in the future.



# [`clean_osm_data_options``

Specifies the options to clean the `OpenStreetMap](https://wiki.osmfoundation.org/wiki/Main_Page) (OSM) data.


   :language: yaml
   :start-at: clean_osm_data_options:
   :end-at: generator_name_method:


   :header-rows: 1
   :widths: 25,10,22,27
   :file: configtables/clean_osm_data_options.csv



# [`build_osm_network``

Specifies the options to build the `OpenStreetMap](https://wiki.osmfoundation.org/wiki/Main_Page) (OSM) network.


   :language: yaml
   :start-at: build_osm_network:
   :end-at: force_ac:


   :header-rows: 1
   :widths: 25,10,22,27
   :file: configtables/build_osm_network.csv



# [`base_network`

Specifies the minimum voltage magnitude in the base network and the offshore substations.


   :language: yaml
   :start-at: base_network:
   :end-at: min_voltage_rebase_voltage:


   :header-rows: 1
   :widths: 25,10,22,27
   :file: configtables/base_network.csv



# `load_options`

Specifies the options to estimate future electricity demand (load). Different years might be considered for weather and the socioeconomic pathway (GDP and population growth), to enhance modelling capabilities.


   :language: yaml
   :start-at: load_options:
   :end-at: scale:


   :header-rows: 1
   :widths: 25,10,22,27
   :file: configtables/load_options.csv


    The snapshots date range (`snapshots\start` - `snapshots\end`) must be in the `weather_year`.



# `co2_budget`

If enabled, this option allows setting different CO₂ targets for each planning horizon year. Only supports foresights with planning horizon such as myopic.


   :language: yaml
   :start-at: co2_budget:
   :end-at: 2050:


   :header-rows: 1
   :widths: 25,10,22,27
   :file: configtables/co2_budget.csv



# `electricity`

Specifies the options for the rule `add_electricity`. This includes options across several features, including but not limited to: voltage levels, electricity carriers available, renewable capacity estimation, CO2 emission limits, operational reserve, storage parameters. See the table below for more details.


   :language: yaml
   :start-at: electricity:
   :end-at: PV:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/electricity.csv


    Carriers in `conventional_carriers` must not also be in `extendable_carriers`.



# `lines`

Specifies electricity line parameters.


   :language: yaml
   :start-after: PV:
   :end-at: under_construction:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/lines.csv



# `links``

Specifies Link parameters. Links are a fundamental component of `PyPSA](https://pypsa.readthedocs.io/en/latest/components.html) .


   :language: yaml
   :start-at: links:
   :end-before: transformers:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/links.csv



# [`transformers`

Specifies transformers parameters and types.


   :language: yaml
   :start-at: transformers:
   :end-before: atlite:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/transformers.csv



# `atlite`

Define and specify the `atlite.Cutout` used for calculating renewable potentials and time-series. All options except for `features`` are directly used as `cutout parameters](https://atlite.readthedocs.io/en/latest/ref_api.html#cutout).


   :language: yaml
   :start-at: atlite:
   :end-before: renewable:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/atlite.csv



# [`renewable`

Specifies the options to obtain renewable potentials in every cutout. These are divided in five different renewable technologies: onshore wind (`onwind`), offshore wind with AC connection (`offwind-ac`), offshore wind with DC connection (`offwind-dc`), solar (`solar`), and hydropower (`hydro`).

## `onwind`


   :language: yaml
   :start-at: renewable:
   :end-before:   offwind-ac:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/onwind.csv

## `offwind-ac`


   :language: yaml
   :start-at:   offwind-ac:
   :end-before:   offwind-dc:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/offwind-ac.csv

## `offwind-dc`


   :language: yaml
   :start-at:   offwind-dc:
   :end-before:   solar:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/offwind-dc.csv

## `solar`


   :language: yaml
   :start-at:   solar:
   :end-before:   hydro:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/solar.csv

## `hydro`


   :language: yaml
   :start-at:   hydro:
   :end-at: multiplier:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/hydro.csv

## `csp`


   :language: yaml
   :start-at:   csp:
   :end-at: csp_model:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/csp.csv



# `costs`

Specifies the cost assumptions of the technologies considered. Cost information is obtained from the config file and the file `data/costs.csv`, which can also be modified manually.


   :language: yaml
   :start-after: Costs Configuration
   :end-at: CCGT: 0.58


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/costs.csv


   To change cost assumptions in more detail (i.e. other than `marginal_cost`), consider modifying cost assumptions directly in `data/costs.csv` as this is not yet supported through the config file.
   You can also build multiple different cost databases. Make a renamed copy of `data/costs.csv` (e.g. `data/costs-optimistic.csv`) and set the variable `COSTS=data/costs-optimistic.csv` in the `Snakefile`.


   The `marginal costs` or in this context `variable costs`` of operating the assets is important for realistic operational model outputs.
   It can define the curtailment order of renewable generators, the dispatch order of generators, and the dispatch of storage units.
   If not approapriate set, the model might output unrealistic results. Learn more about this in
   `Parzen et al. 2023](https://www.sciencedirect.com/science/article/pii/S2589004222020028) and in
   [Kittel et al. 2022](https://www.sciencedirect.com/science/article/pii/S2589004222002723).


# `monte_carlo`

Specifies the options for Monte Carlo sampling.


   :language: yaml
   :start-at: monte_carlo:
   :end-before:   solving:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/monte-carlo.csv



# `solving`

Specify linear power flow formulation and optimization solver settings.

## `options`


   :language: yaml
   :start-at: solving:
   :end-before:   solver:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/solving-options.csv

## `solver`


   :language: yaml
   :start-at:   solver:
   :end-before: plotting:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/solving-solver.csv



# `plotting`

Specifies plotting options.


   :language: yaml
   :start-at: plotting:


   :header-rows: 1
   :widths: 25,7,22,30
   :file: configtables/plotting.csv
