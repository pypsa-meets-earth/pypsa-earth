<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

# Configuration

PyPSA-Earth imports the configuration options originally developed in [PyPSA-Eur](https://pypsa-eur.readthedocs.io/en/latest/index.html) and here reported and adapted.
The options here described are collected in a `config.yaml` file located in the root directory.
Users should copy the provided default configuration (`config.default.yaml`) and amend
their own modifications and assumptions in the user-specific configuration file (`config.yaml`);
confer installation instructions at [installation](../home/installation.md).

  Credits to PyPSA-Eur developers for the initial drafting of the configuration documentation here reported

## Keeping your config up to date

PyPSA-Earth loads ``config.default.yaml`` first and then merges your ``config.yaml`` on top (see the ``configfile`` entries in the Snakefile). You therefore only need to list settings in ``config.yaml`` that differ from the defaults—a small override file is enough for most studies.

When you upgrade to a new version, open the updated ``config.default.yaml`` and check what changed: new keys, renamed paths, or reorganised sections. Compare it with your ``config.yaml`` and copy across any new defaults you want to use, or move keys you still override to their new locations. The [release notes](../release-notes.md) summarise breaking config changes per release.

Some deprecated key names are migrated automatically when the workflow starts; if Snakemake prints a ``FutureWarning`` about an old path, update your ``config.yaml`` to the new key and remove the obsolete one. See [Renamed keys](#renamed-keys) for the full mapping. For study-specific settings you can also pass extra files with ``snakemake --configfile my_study.yaml`` instead of growing a single ``config.yaml``.

### Renamed keys

The table below lists all keys that have been renamed or moved. The old keys still work (their values are copied to the new location), but you should update your ``config.yaml`` to avoid future breakage.

| Old key (deprecated) | New key |
|---|---|
| `electricity.co2limit` | `co2.limit` |
| `electricity.co2base` | `co2.base` |
| `electricity.automatic_emission` | `co2.automatic_emission.enable` |
| `electricity.automatic_emission_base_year` | `co2.automatic_emission.base_year` |
| `costs.emission_prices.co2` | `co2.emission_price` |
| `co2_budget.enable` | `co2.budget.enable` |
| `co2_budget.override_co2opt` | `co2.budget.override_co2opt` |
| `co2_budget.year` | `co2.budget.year` |
| `co2_budget.co2base_value` | `co2.budget.base_value` |
| `sector.solar_thermal` *(bool flag)* | `sector.solar_thermal_collector.enable` |
| `sector.solar_cf_correction` | `sector.solar_thermal_collector.cf_correction` |
| `solar_thermal.clearsky_model` | `sector.solar_thermal_collector.clearsky_model` |
| `solar_thermal.orientation` | `sector.solar_thermal_collector.orientation` |
| `clean_osm_data_options` | `osm.clean_osm_data` |
| `build_osm_network` | `osm.build_osm_network` |
| `fossil_reserves.{carrier}` | `sector.{carrier}.reserves` |

`{carrier}` is the fuel name (e.g. `oil`, `coal`, `gas`, `lignite`, `biomass`). Migrations run automatically via ``migrate_config`` in ``scripts/_helpers.py``; see also the [release notes](../release-notes.md) when upgrading.

## Top-level configuration

Version, tutorial mode, logging, and optional custom Snakemake rules (``# META`` banner in ``config.default.yaml``).

```yaml
--8<-- "configtables/snippets/meta.yaml"
```

{{ read_csv('configtables/meta.csv') }}

## Study setup

Options under the **STUDY SETUP** banner in ``config.default.yaml``: region, run namespacing, [scenario wildcards](wildcards.md), model time range, and output directories.

Region, foresight mode, and default output directories (`results_dir`, `summary_dir` at the end of the study setup block, not nested under ``run``).

```yaml
--8<-- "configtables/snippets/study_setup.yaml"
```

{{ read_csv('configtables/study_setup.csv') }}

### run

It is common conduct to analyse energy system
optimisation models for **multiple scenarios** for
a variety of reasons,
e.g. assessing their sensitivity towards changing
the temporal and/or geographical resolution or
investigating how
investment changes as more ambitious greenhouse-gas
emission reduction targets are applied.

The `run` section is used for running and storing
scenarios with different configurations which are
not covered by [wildcards](wildcards.md). It
determines the path at which resources, networks
and results are stored. Therefore the user can run
different configurations within the same directory.
If a run with a non-empty name should use cutouts
shared across runs, set `shared_cutouts` to `true`.

```yaml
--8<-- "configtables/snippets/run.yaml"
```

{{ read_csv('configtables/run.csv') }}

### scenario

The `scenario` section is an extraordinary section of the config file
that is strongly connected to the [wildcards](wildcards.md) and is designed to
facilitate running multiple scenarios through a single command

```bash
snakemake -j 1 solve_all_networks
```

For each wildcard, a **list of values** is provided. The rule `solve_all_networks` will trigger the rules for creating `results/networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc` for **all combinations** of the provided wildcard values as defined by Python's [itertools.product(...)](https://docs.python.org/2/library/itertools.html#itertools.product) function that snakemake's [expand(...) function](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#targets) uses.

An exemplary dependency graph (starting from the simplification rules) then looks like this:

![Image](https://raw.githubusercontent.com/pypsa-meets-earth/documentation/main/doc/img/scenarios.png)

```yaml
--8<-- "configtables/snippets/scenario.yaml"
```

{{ read_csv('configtables/scenario.csv') }}

### snapshots

Specifies the temporal range for the historical weather data, which is used to build the energy system model. It uses arguments to [pandas.date_range](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.date_range.html). The date range must be in the past (before 2022). A well-tested year is 2013.

```yaml
--8<-- "configtables/snippets/snapshots.yaml"
```

{{ read_csv('configtables/snapshots.csv') }}

## Data retrieval

Switches in ``enable`` for downloading databundles, OSM data, cutouts, and building inputs locally. Set most retrieve flags to ``true`` for the first run.

```yaml
--8<-- "configtables/snippets/data_retrieval.yaml"
```

{{ read_csv('configtables/data_retrieval.csv') }}

## Geography & shapes

Coordinate systems, regional shapes, subregions, land-cover exclusions, OpenStreetMap, clustering, and optional augmented connectivity (``# GEOGRAPHY & SHAPES`` in ``config.default.yaml``).

### crs

Defines the coordinate reference systems (crs).

```yaml
--8<-- "configtables/snippets/crs.yaml"
```

{{ read_csv('configtables/crs.csv') }}

### build_shape_options

Specifies the options to build the shapes in
which the region of interest (`countries`) is
divided.

```yaml
--8<-- "configtables/snippets/build_shape_options.yaml"
```

{{ read_csv('configtables/build_shape_options.csv') }}

### subregion

If enabled, this option allows a region of interest (`countries`) to be redefined into subregions,
which can be activated at various stages of the workflow. Currently, it is used in `simplify_network` and `cluster_network` rule.

```yaml
--8<-- "configtables/snippets/subregion.yaml"
```

{{ read_csv('configtables/subregion.csv') }}

The names of subregions are arbitrary. Its sizes are determined by how many GADM IDs that are included in the list.
A single country can be divided into multiple subregions, and a single subregion can include GADM IDs from multiple countries.
If the same GADM ID appears in different subregions, the first subregion listed will take precedence over that region.
The remaining GADM IDs that are not listed will be merged back to form the remaining parts of their respective countries.
For example, consider the Central District of Botswana, which has a GADM ID of `BW.3`. To separate this district from the rest of the country, you can select:

> **See `config.default.yaml` for the full configuration.**

There are several formats for GADM IDs depending on the version, so before using this feature, please review the `resources/shapes/gadm_shape.geojson` file which can be created using the command:

```bash
snakemake -j 1 build_shapes
```

!!! note
    The rule `build_shapes` currently use [Version 4.1](https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/) for their GADM data. This may change in the future.

### natura

Land-cover exclusion raster settings. Only relevant when ``enable.build_natura_raster`` is ``true``.

```yaml
--8<-- "configtables/snippets/natura.yaml"
```

{{ read_csv('configtables/natura.csv') }}

### osm

OpenStreetMap ([OSM](https://wiki.osmfoundation.org/wiki/Main_Page)) pipeline settings live under the ``osm`` block: ``osm.clean_osm_data`` filters raw downloads; ``osm.build_osm_network`` merges buses and lines into the base-network CSV inputs. If your config still uses the old top-level names, see [Renamed keys](#renamed-keys).

```yaml
--8<-- "configtables/snippets/osm.yaml"
```

{{ read_csv('configtables/osm.csv') }}

### cluster_options

Specifies the options to simplify and cluster the
network. This is done in two stages, first using
the rule `simplify_network` and then using the
rule `cluster_network`. For more details on this
process, see the [PyPSA-Earth paper](https://www.
sciencedirect.com/science/article/pii/
S0306261923004609), section 3.7.

```yaml
--8<-- "configtables/snippets/cluster_options.yaml"
```

{{ read_csv('configtables/cluster_options.csv') }}

### augmented_line_connection

If enabled, it increases the connectivity of the
network. It makes the network graph
[k-edge-connected](https://en.wikipedia.org/wiki/
K-edge-connected_graph), i.e.,
if fewer than k edges are removed, the network
graph stays connected. It uses the
[k-edge-augmentation](https://networkx.org/
documentation/stable/reference/algorithms/
generated/networkx.algorithms.connectivity.
edge_augmentation.k_edge_augmentation.
html#networkx.algorithms.connectivity.
edge_augmentation.k_edge_augmentation)
algorithm from the [NetworkX](https://networkx.
org/documentation/stable/index.html) Python
package.

```yaml
--8<-- "configtables/snippets/augmented_line_connection.yaml"
```

{{ read_csv('configtables/augmented_line_connection.csv') }}


## Network & resources

Base network assembly, demand, electricity grid, weather and renewables, costs and emissions, and Monte Carlo settings (``# NETWORK & RESOURCES`` in ``config.default.yaml``).

### base_network

Specifies the minimum voltage magnitude in the base network and the offshore substations.

```yaml
--8<-- "configtables/snippets/base_network.yaml"
```

{{ read_csv('configtables/base_network.csv') }}

### load_options

Specifies the options to estimate future electricity demand (load). Different years might be considered for weather and the socioeconomic pathway (GDP and population growth), to enhance modelling capabilities.

```yaml
--8<-- "configtables/snippets/load_options.yaml"
```

{{ read_csv('configtables/load_options.csv') }}

The snapshots date range (`snapshots\start` - `snapshots\end`) must be in the `weather_year`.

### electricity

Specifies the options for the rule `add_electricity`. This includes options across several features, including but not limited to: voltage levels, electricity carriers available, renewable capacity estimation, operational reserve, storage parameters. See the table below for more details.

```yaml
--8<-- "configtables/snippets/electricity.yaml"
```

{{ read_csv('configtables/electricity.csv') }}

Carriers in `conventional_carriers` must not also be in `extendable_carriers`.

### lines

Specifies electricity line parameters.

```yaml
--8<-- "configtables/snippets/lines.yaml"
```

{{ read_csv('configtables/lines.csv') }}

### links

Specifies Link parameters. Links are a fundamental component of [PyPSA](https://pypsa.readthedocs.io/en/latest/components.html) .

```yaml
--8<-- "configtables/snippets/links.yaml"
```

{{ read_csv('configtables/links.csv') }}

### transformers

Specifies transformers parameters and types.

```yaml
--8<-- "configtables/snippets/transformers.yaml"
```

{{ read_csv('configtables/transformers.csv') }}

### atlite

Define and specify the `atlite.Cutout` used for calculating renewable potentials and time-series. All options except for `features` are directly used as [cutout parameters](https://atlite.readthedocs.io/en/latest/ref_api.html#cutout).

{{ read_csv('configtables/atlite.csv') }}

### renewable

Specifies the options to obtain renewable potentials in every cutout. These are divided in five different renewable technologies: onshore wind (`onwind`), offshore wind with AC connection (`offwind-ac`), offshore wind with DC connection (`offwind-dc`), solar (`solar`), and hydropower (`hydro`).

#### onwind

```yaml
--8<-- "configtables/snippets/renewable_onwind.yaml"
```

{{ read_csv('configtables/onwind.csv') }}

#### offwind-ac

```yaml
--8<-- "configtables/snippets/renewable_offwind-ac.yaml"
```

{{ read_csv('configtables/offwind-ac.csv') }}

#### offwind-dc

```yaml
--8<-- "configtables/snippets/renewable_offwind-dc.yaml"
```

{{ read_csv('configtables/offwind-dc.csv') }}

#### solar

```yaml
--8<-- "configtables/snippets/renewable_solar.yaml"
```

{{ read_csv('configtables/solar.csv') }}

#### hydro

```yaml
--8<-- "configtables/snippets/renewable_hydro.yaml"
```

> **See `config.default.yaml` for the full configuration.**

{{ read_csv('configtables/hydro.csv') }}

#### csp

```yaml
--8<-- "configtables/snippets/renewable_csp.yaml"
```

> **See `config.default.yaml` for the full configuration.**

{{ read_csv('configtables/csp.csv') }}

### costs

Specifies the cost assumptions of the technologies considered. Cost information is obtained from the config file and the file `data/costs.csv`, which can also be modified manually.

```yaml
--8<-- "configtables/snippets/costs.yaml"
```

{{ read_csv('configtables/costs.csv') }}

!!! note
    To change cost assumptions in more detail (i.e. other than `marginal_cost`), consider modifying cost assumptions directly in `data/costs.csv` as this is not yet supported through the config file.
    You can also build multiple different cost databases. Make a renamed copy of `data/costs.csv` (e.g. `data/costs-optimistic.csv`) and set the variable `COSTS=data/costs-optimistic.csv` in the `Snakefile`.

    The `marginal costs` or in this context `variable costs` of operating the assets is important for realistic operational model outputs.
    It can define the curtailment order of renewable generators, the dispatch order of generators, and the dispatch of storage units.
    If not appropriately set, the model might output unrealistic results. Learn more about this in
    [Parzen et al. 2023](https://www.sciencedirect.com/science/article/pii/S2589004222020028) and in
    [Kittel et al. 2022](https://www.sciencedirect.com/science/article/pii/S2589004222002723).

### co2

Carbon dioxide policy settings: emission caps and carbon prices (via `Co2L`/`Ep` wildcards), automatic emission extraction, and optional planning-horizon budgets in sector-coupled runs.

```yaml
--8<-- "configtables/snippets/co2.yaml"
```

{{ read_csv('configtables/co2.csv') }}

### monte_carlo

Specifies the options for Monte Carlo sampling.

```yaml
--8<-- "configtables/snippets/monte_carlo.yaml"
```

{{ read_csv('configtables/monte-carlo.csv') }}

## Sector options

Options under the **SECTOR OPTIONS** banner in ``config.default.yaml``: hydrogen policy and export, sector demand inputs, custom data overrides, brownfield capacity grouping, and sector coupling settings.

### policy_config

Specifies energy-policy options for hydrogen, for example temporal matching and additionality constraints.

```yaml
--8<-- "configtables/snippets/policy_config.yaml"
```

{{ read_csv('configtables/policy_config.csv') }}

### export

Specifies hydrogen export demand, storage, and shipping profile settings.

```yaml
--8<-- "configtables/snippets/export.yaml"
```

{{ read_csv('configtables/export.csv') }}

### demand_data

Specifies sector-coupled demand inputs (UNSD energy balances and related years).

```yaml
--8<-- "configtables/snippets/demand_data.yaml"
```

{{ read_csv('configtables/demand_data.csv') }}

### custom_data

Specifies which custom datasets are used to replace or supplement the default model data. For full details see [Custom Data Integration](custom-data.md).

```yaml
--8<-- "configtables/snippets/custom_data.yaml"
```

{{ read_csv('configtables/custom_data.csv') }}

### existing_capacities

Vintage grouping and thresholds for brownfield power and heating capacities in ``add_existing_baseyear`` and related rules.

```yaml
--8<-- "configtables/snippets/existing_capacities.yaml"
```

{{ read_csv('configtables/existing_capacities.csv') }}

### sector

Specifies the options for sector coupling, i.e. the integration of the electricity system with other sectors such as heating and transport.

#### top-level

Carrier toggles, fossil-fuel supply settings (`gas`, `coal`, `lignite`, `oil`), hydrogen, and ammonia. Fossil fuel reserves are set per carrier as `sector.{carrier}.reserves` [TWh/bus] (e.g. `sector.oil.reserves`, `sector.coal.reserves`). The value sets initial Store energy in `add_carrier_buses` for fuel carriers used in `sector.conventional_generation` (`gas`, `oil`, `coal`, `lignite`, `biomass`); it defaults to 0 if omitted. The former top-level ``fossil_reserves`` block is deprecated — see [Renamed keys](#renamed-keys).

```yaml
--8<-- "configtables/snippets/sector_toplevel.yaml"
```

{{ read_csv('configtables/sector_toplevel.csv') }}

#### heat sector

Solar thermal collector settings live under ``sector.solar_thermal_collector`` (not ``sector.solar_thermal``). Deprecated solar-thermal keys are listed in [Renamed keys](#renamed-keys).

```yaml
--8<-- "configtables/snippets/sector_heat.yaml"
```

{{ read_csv('configtables/sector_heat.csv') }}

#### land transport sector

```yaml
--8<-- "configtables/snippets/sector_land_transport.yaml"
```

{{ read_csv('configtables/sector_land_transport.csv') }}

#### biomass sector

```yaml
--8<-- "configtables/snippets/sector_biomass.yaml"
```

{{ read_csv('configtables/sector_biomass.csv') }}

#### electricity distribution grid

```yaml
--8<-- "configtables/snippets/sector_electricity_distribution_grid.yaml"
```

{{ read_csv('configtables/sector_electricity_distribution_grid.csv') }}

#### shipping & aviation sector

```yaml
--8<-- "configtables/snippets/sector_shipping_aviation.yaml"
```

{{ read_csv('configtables/sector_shipping_aviation.csv') }}

#### ccus & conversion options

```yaml
--8<-- "configtables/snippets/sector_ccus.yaml"
```

{{ read_csv('configtables/sector_ccus.csv') }}

#### industry options

```yaml
--8<-- "configtables/snippets/sector_industry.yaml"
```

{{ read_csv('configtables/sector_industry.csv') }}

#### powerplants options

```yaml
--8<-- "configtables/snippets/sector_powerplants.yaml"
```

{{ read_csv('configtables/sector_powerplants.csv') }}

## Solving

Options under the **SOLVING** banner in ``config.default.yaml``: solver choice, linear formulation, load shedding, iteration settings, solver presets, and memory limit.

### solver

```yaml
--8<-- "configtables/snippets/solving_solver.yaml"
```

{{ read_csv('configtables/solving-solver.csv') }}

### options

```yaml
--8<-- "configtables/snippets/solving_options.yaml"
```

{{ read_csv('configtables/solving-options.csv') }}

## Plotting

Options under the **PLOTTING** banner in ``config.default.yaml``: map layout, plot thresholds, technology groupings, carrier colours, and display names.

```yaml
--8<-- "configtables/snippets/plotting.yaml"
```

{{ read_csv('configtables/plotting.csv') }}
