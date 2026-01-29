<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

Release Notes

# Upcoming release

This part of documentation collects descriptive release notes to capture the main improvements introduced by developing the model before the next release.

**New Features and Major Changes**

* Drop use of override_components that is no longer needed in newer PyPSA versions [PR #1699](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1699)

* Add an option to redefine countries into subregions in `cluster_networks` [PR #1542](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1542)

* Added Section covering Optimization for Energy Systems Models [PR #1558](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1558)

* Add the configuration `co2_budget` to set CO₂ emission targets in multiple planning horizon years [PR #1553](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1553)

* Add pseudo `branch()` to streamline snakemake workflow. Replace this with Snakemake implementation if Snakemake >= 8.3.0. Added an option to configure which sector components to include which when disabled, irrelevant rules are skipped automatically [PR #1538](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1538)

* Extend functionality of isolated nodes to networks and fix to demand scaling when gdp or pop is empty [PR #1540](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1540)

* Redesign temporal matching including yearly, monthly and hourly matching [PR #1463](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1463)

* Add option to ignore loading network data in clean_osm_data [PR #1580](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1580)

* Add ll wildcard option 'l{factor}' to enable line-wise transmission expansion, including optional lower and upper bounds [PR #1592](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1592)

**Minor Changes and bug-fixing**

* Bump powerplantmatching to 0.8.0 [PR #1702](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1702)

* Update reference values of the objective function in validator workflow and adjust format of the objective outputs in csv files [PR #1705](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1705)

* Update steel GEM data to version 2024 and use backup link for GEM pipelines [PR #1708](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1708)

* Refine load shedding capacity calculation to use bus-specific maximum loads instead of fixed large values, improving solver performance and numeric stability. Rename load shedding generator carrier to 'load shedding' (now) instead of just 'load' (previously) [PR #1581](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1581)

* The configuration setting for ``focus_weights` has been moved from [focus_weights:` to `cluster_options: focus_weights:``. Backwards compatibility to old config files is maintained `PR #1565](<https://github.com/pypsa-meets-earth/pypsa-earth/pull/1565>)

* Disable distribute_cluster only for one subnetwork and country [PR #1539](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1539)

* Reintroduce sanitize_carriers and sanitize_location to reduce the number of warnings related to components with undefined carriers [PR #1555](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1555)

* Fix the offwind depth calculation by providing the proper GEBCO file [PR #1559](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1559)

* bug-fix in `clean_osm_data.py` so that cleaned data has the same number of circuits as in the raw data and assumptions are correctly applied [PR #1552](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1552)

* Avoid using variable "nodes" in prepare_sector_network unless explicitly in the arguments [PR #1575](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1575)

* Fix and reactivate the option for a custom busmap in `cluster_network.py` [PR #1537](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1537)

* Bug-fixing override co2opt in add_co2_budget [PR #1597](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1597)

* Fix shipping and aviation implementation in multi-country models [PR #1582](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1582)

* Add missing colors for energy carriers in the sector-coupled model [PR #1625](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1625)

* Fix Labeling of technologies [PR #1644](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1644)

* fix workflow when augmmented_line_connection is false for Tunisia [PR #1677](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1677)

# PyPSA-Earth 0.7.0

**New Features and Major Changes**

* Revise bundle_cutouts_northamerica [PR #1479](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1479)

* Revise pinned to locked environments [PR #1458](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1458)

* Update to `technology-data` version v0.12.0 [PR #1452](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1452)

* Add a command-line interface in `databundle_cli.py` that will be triggered if `retrieve_databundle_light.py` fails to retrieve all necessary files, providing a fallback to assist with debugging the issue [PR #1366](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1366)

* Add Wikipedia as a source for the preparation of transport_data.csv  [PR #1410](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1410)

* Add a new page for sector-coupled tutorial [PR #1374](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1374)

* Add an option to redefine countries into subregions which can be activated at various stages of the workflow. The subregions can be defined either by the GADM_ID or a custom map. Currently, it is only used in the `simplify_network` rule [PR #1300](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1300)

* Drop duplication of retrieve_data and COST_DIR, add params and update technology-data version [PR #1249](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1249)

* In alternative clustering, generate hydro inflows by shape and avoid hydro inflows duplication for plants installed in the same node [PR #1120](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1120)

* Add a function to calculate length-based efficiencies and apply it to the H2 pipelines [PR #1192](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1192)

* Support of Linopy for Power and Sector-Coupled Modelling and latest PyPSA version [PR #1172](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1172)

* Update workflow to geopandas >= 1.0 [PR #1276](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1276)

* Index hydro units by their location and drop use of alternative clustering option for hydro [PR #1331](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1331)

* Introduce universal currency conversion to allow use of currencies other than EUR in input data and results [PR #1319](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1319)

* Introduce US-specific cost configurations and update to `technology_data_version` v0.13.2 [PR #1448](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1448)

* Add new hydrogen production technologies and reorganize existing structure [PR #1227](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1227)

**Minor Changes and bug-fixing**

* Remove duplicates in environment file [PR #1473](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1473)

* Fix fallback logic for missing currency data by correctly applying the default exchange rate and build cache for storing commonly used exchange rates [PR #1492](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1492)

* Update code to reflect the Wikipedia source page for `prepare_transport_data_input` [PR #1486](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1486)

* Add efficiency gain and growth rates for other energy consumption and fill missing NaNs with 0 [PR #1468](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1468)

* Fix missing nodes in prepare_sector_network [PR #1432](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1432)

* Fix missing country attribute in cluster_network [PR #1443](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1443)

* Fix numpy TypeErrors [PR #1430](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1430)

* Update description of the sector-coupled capabilities [PR #1449](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1449)

* Modified how_to_contribute documentation page [PR #1439](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1439)

* Redefine h2export wildcard_constraint to allow digits + optional dot followed by more digites (decimals) [PR #1434](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1434)

* Update the AL_production.csv data from another source with more countries [PR #1428](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1428)

* Fix params for prepare_sector_network script [PR #1427](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1427)

* Add option to adjust load shedding costs [PR #1403](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1403)

* Fix: Use the mean value instead of the sum to remove duplicates in the urban percentage data per country and year [PR #1420](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1420)

* Fix prepare_ports script to generate both outputs if export ports with custom data is selected [PR #1424](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1424)

* Update the applications list for PyPSA-Earth model [PR #1413](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1413)

* Fix problem with a discontinued World Bank data source in prepare_transport_data_input [PR #1410](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1410)

* Fix bug in myopic run [PR #1369](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1369)

* Fix missing focus_weights on cluster_network params without augmented_line_connections, minimize warnings on subregions [PR #1348](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1348)

* Align structure of the components with consistency checks in the updated PyPSA version [PR #1315](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1315)

* Prevent computation of powerplantmatching if replace option is selected for custom_powerplants [PR #1281](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1281)

* Fix overlapping bus regions when alternative clustering is selected [PR #1287](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1287)

* Fix readthedocs by explicitly specifying the location of the Sphinx config [PR #1292](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1292)

* Fix lossy bidirectional links, especially H2 pipelines, which would sometimes gain H2 instead of losing it [PR #1192](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1192)

* Fix the need for administrative rights on Windows by changing all shadow directory settings for Windows in the Snakefile [PR #1295](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1295) and  [PR #1301](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1301)

* Add CITATION.cff to guide users on how cite PyPSA-Earth. Insert the DOI badge in the README linking to the very first PyPSA-Earth paper [PR #1316](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1316)

* Remove pyomo from cluster_network and introduce scip [PR #1320](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1323)

* Fix mixup of plant locations for hydro inflow and inflow overestimation [PR #1322](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1322)

* Fix namibia geofk, line country tag mismatch and minor fixes [PR #1330](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1330)

* Minor revision docker workflow to have it working on upstream [PR #1343](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1343)

* Fix the scaling factor for time-varying loads of the sector model [PR #1372](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1372)

* Integrate EIA data for US-specific CAGR and fuel shares for the sector model [PR #1372](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1400)

* Revise naming of Wikipedia data for vehicles [PR #1422](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1422)

* Monte Carlo: move qmc.discrepancy to logging only [PR #1418](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1418)

* Extending powerplant filter option to custom powerplants [PR #1465](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1465)

* Fix an issue with the GEBCO file by limiting libgdal-core<3.10 [PR #1519](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1519)

* Fix a naming issue with European cutout [PR #1530](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1530)

* Avoid adding CO2 pipeline links when option is disabled [PR #1504](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1504)

* Add possibility to overwrite cost attributes for sector model [PR #1567](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1567)

# PyPSA-Earth 0.6.0

**New Features and Major Changes (24th December 2024)**

* Include option in the config to allow for custom airport data [PR #1241](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1241)

* Added Dev Containers and docker as an option to get started with pypsa-earth [PR #1228](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1228)

* Add a list of PyPSA-Earth applications in academic and industrial projects [PR #1255](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1255)

* Computational improvements of build_osm_network [PR #845](https://github.com/pypsa-meets-earth/pypsa-earth/pull/845)

* Boost computational performances of set_lines_ids with cKDTree by scipy [PR #806](https://github.com/pypsa-meets-earth/pypsa-earth/pull/806)

* Boost computational performances of set_substation_ids using DBSCAN [PR #799](https://github.com/pypsa-meets-earth/pypsa-earth/pull/799)

* Boost computational performances of fix_overpassing_line [PR #807](https://github.com/pypsa-meets-earth/pypsa-earth/pull/807)

**Minor Changes and bug-fixing**

* Added electricity bus to Fischer-Tropsch in prepare_sector_network.py [PR #1226](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1226)

* Update BW, NG and BJ tutorial databundles to include gadm-like sources from geoboundaries [PR #1257](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1257)

# PyPSA-Earth 0.5.0

**New Features and Major Changes (14th December 2024)**

* Added capabilities of cross-sectoral modeling by merging with PyPSA-Earth-Sec model [https://github.com/pypsa-meets-earth/pypsa-earth-sec`__

* The workflow configuration now supports incremental changes to the default configuration in the `config.yaml` and configfiles passed to snakemake via `--configfile myconfig.yaml`. Therefore the user may now only include settings in their `config.yaml` which differ from the default configuration. One can think of the new `config.yaml` as of a list of arguments in a python function that already have a default. So in principle the `config.yaml` could now be empty, and the workflow would still run. `PR #1053](<https://github.com/pypsa-meets-earth/pypsa-earth/pull/1053>)

* Include option of endogenous export, which optimizes the export quantity based on price signals [PR #1201](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1201)

* Remove elec-based H2 and battery technologies before addition in `prepare_sector_network.py` script and fix bus names for links that models H2 repuspose network [PR #1198](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1198)

* Add electricity distribution grid with solar rooftop and home battery technologies [PR #1221](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1221)

* Include a dedicated cutout for North America in bundle_config.yaml [PR #1121](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1121)

* Include a dedicated cutout for Europe in bundle_config.yaml [PR #1125](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1125)

* Include a dedicated cutout for Oceania in bundle_config.yaml [PR #1157](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1157)

* Integrate RDIR into sector rules to store intermediate data in scenario folders [PR #1154](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1154)

* The computation of `hydro_profile.nc` in [build_renewable_profiles.py` is not differentiated whether alternative clustering is applied or not; the indexing of the different power plants in `add_electricity.py` is performed according to the bus either in case alternative clustering is applied or not and a `hydro_inflow_factor` is computed prior to the computation of `inflow_t` to split the inflow according to the capacity of each different unit of each power plant (if more units are present). `PR #1119](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1119)

* Use BASE_DIR in rules and `_helpers.py` script for facilitate module import in subworkflow [PR #1137](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1137)

* Enable sector rules import in subworkflow [PR #1178](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1178)

**Minor Changes and bug-fixing**

* The default configuration for `electricity:estimate_renewable_capacities:year` was updated from 2020 to 2023. [PR #1106](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1106)

* Fix the mismatch between buses and x, y locations while creating H2 Stores [PR #1134](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1134)

* Enable configfile specification for mock_snakemake [PR #1135](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1135)

* Removed duplications of devendencies in environment.yaml [PR #1128](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1128)

* Fix pre-commit docformatter python issue. [PR #1153](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1153)

* Drop duplicate entries in `AL_production.csv` data used in [build_industry_demand` rule `PR #1143](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1143)

* Fix bugs in `prepare_sector_network.py` related to links with H2 buses and bug of re-addition of H2 and battery carriers in present [PR #1145](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1145)

* Drop entries that contain non-string elements in country column of `CO2_emissions_csv` data in [prepare_transport_data_input.py` script `PR #1166](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1166)

* Local tests are now run with `make test`. This uses a `Makefile` which runs snakemake calls with different configurations. [PR #1053](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1053)

* Adds [Dependabot](https://docs.github.com/en/code-security/getting-started/dependabot-quickstart-guide) to keep GitHub actions up to date. [PR #1184](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1184)

* Adds code security scans via [CodeQL](https://codeql.github.com/) to CI. [PR #1185](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1185)

* Adds CI to update keep pinned environment files up to date. [PR #1183](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1183) and  [PR #1210](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1210)

* Revise ports data for export in `add_export.py` related to sector model [PR #1175](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1175)

* Restore string values of tech_colors in config file [PR #1205](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1205)

* Drop vrestil dependency [PR #1220](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1220)

* Include a configuration option to simplify / not simplify shapefiles based on a boolean value specified under [build_shape_options:simplify_gadm` option in the config file `PR 1138](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1138)

* Fix the mismatch between buses and x, y locations while creating H2 Stores [PR #1134](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1134)

* Remove duplicate entries from hydrogen export ports [PR #1233](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1233)

* Fix the environment placing a version limit to numpoly [PR #1237](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1237)

# PyPSA-Earth 0.4.1

**New Features and Major Changes (19th September 2024)**

* Add functionality to modify the cost assumptions using config files [PR #1097](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1097)

**Minor Changes and bug-fixing**

* Remove unused `countries_codes` argument from [load_GDP` function in `build_shapes.py` script, which was not being called as intended with positional arguments `PR #1069](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1069)
* Fixed problematic float parsing (`_parse_float`) in `clean_osm_data.py` to make sure all OSM lines are correctly accounted for [PR #1089](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1089)
* Fix minor bug for advanced csp implementation [PR #1076](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1076)
* Fix minor bug in `build_powerplants.py` where the gas technology assignment incorrectly introduced NaN values for all powerplant technologies. [PR #1102](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1102)

# PyPSA-Earth 0.4.0

**New Features and Major Changes (27th July 2024)**

* Improve Monte Carlo feature with more distributions types, independent by PyPSA component. [PR #930](https://github.com/pypsa-meets-earth/pypsa-earth/pull/930)

* Introduce flexible regional selection of the demand files of GEGIS. [PR #991](https://github.com/pypsa-meets-earth/pypsa-earth/pull/991)

* Generalize line types for AC and DC networks. [PR #999](https://github.com/pypsa-meets-earth/pypsa-earth/pull/999)

* Add an option to merge isolated networks into respective backbone networks by countries. [PR #903](https://github.com/pypsa-meets-earth/pypsa-earth/pull/903)

* Add an option to use csv format for custom demand imports. [PR #995](https://github.com/pypsa-meets-earth/pypsa-earth/pull/995)

**Minor Changes and bug-fixing**

* Minor bug-fixing to run the cluster wildcard min [PR #1019](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1019)

* Add option to adjust load scale for each individual countries [PR #1006](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1006)

* Minor bug-fixing to get the generalised line types work for DC lines and AC lines. [PR #1008](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1008) , [PR #1011](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1011) and [PR #1013](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1013)

* Minor bug-fixing for GADM_ID format naming. [PR #980](https://github.com/pypsa-meets-earth/pypsa-earth/pull/980), [PR #986](https://github.com/pypsa-meets-earth/pypsa-earth/pull/986) and [PR #989](https://github.com/pypsa-meets-earth/pypsa-earth/pull/989)

* Fix download_osm_data compatibility for earth-osm v2.1. [PR #954](https://github.com/pypsa-meets-earth/pypsa-earth/pull/954) and [PR #988](https://github.com/pypsa-meets-earth/pypsa-earth/pull/988)

* Improve geometry filtering in clean_osm_data. [PR #989](https://github.com/pypsa-meets-earth/pypsa-earth/pull/989)

* Revise bus region definition by gadm. [PR #1001](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1001)

* Documentation improvements. [PR #1007](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1007)

* Remove unnecessary imports. [PR #1020](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1020)

* Resolve pandas deprecation warning. [PR #1023](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1023)

* Create files where the code outputs the value of the objective function. [PR #1033](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1033)

* Introduce versioning of the configuration files. [PR #1058](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1058)

* Fix bug for hydro inflow normalization for gadm regions (alternative clustering). [PR #1057](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1057)

* Minor bug-fixing for s_nom_min. [PR #961](https://github.com/pypsa-meets-earth/pypsa-earth/pull/961)

# PyPSA-Earth 0.3.0

**New Features and major Changes (24th December 2023)**

* Keep all traceback in logs. [PR #898](https://github.com/pypsa-meets-earth/pypsa-earth/pull/898)

* Function added in clean_osm_data script to allow the use of custom network data instead or on-top of OSM data. [PR #842](https://github.com/pypsa-meets-earth/pypsa-earth/pull/842)

* Improve retrieve_databundle to prioritize smallest databundles [PR #911](https://github.com/pypsa-meets-earth/pypsa-earth/pull/911)

* Add functionality to load shapefiles for hydrobasins directly from the data source directly [PR #919](https://github.com/pypsa-meets-earth/pypsa-earth/pull/919)

* Use [new CC0 v1 dataset](https://doi.org/10.7910/DVN/XIV9BL) for the natura input and automate download of WDPA protected planet data [PR #913](https://github.com/pypsa-meets-earth/pypsa-earth/pull/913)

**Minor Changes and bug-fixing**

* Revise databundles and improve logging in retrieve_databundle [PR #928](https://github.com/pypsa-meets-earth/pypsa-earth/pull/928)

* Improve documentation on installation and short tutorial [PR #918](https://github.com/pypsa-meets-earth/pypsa-earth/pull/918)

# PyPSA-Earth 0.2.3

**New Features and major Changes (19th October 2023)**

* Add params: section in rule definition to keep track of changed settings in config.yaml. [PR #823](https://github.com/pypsa-meets-earth/pypsa-earth/pull/823) and [PR #880](https://github.com/pypsa-meets-earth/pypsa-earth/pull/880)

* Fix Natural Gas implementation in "add_electricity" to avoid "Natural Gas" to be filtered out [PR #797](https://github.com/pypsa-meets-earth/pypsa-earth/pull/797)

* Improve network simplification routine to account for representation HVDC as Line component [PR #743](https://github.com/pypsa-meets-earth/pypsa-earth/pull/743)

* Remove deprecated pypsa.networkclustering approach and replace by pypsa.clustering.spatial [PR #786](https://github.com/pypsa-meets-earth/pypsa-earth/pull/786)

* Drop code-dependency from vresutil [PR #803](https://github.com/pypsa-meets-earth/pypsa-earth/pull/803)

* Add a check to ensure match between a cutout and a modelled area [PR #805](https://github.com/pypsa-meets-earth/pypsa-earth/pull/805)

* Support renewables or renewable expansion to meet a desired share of total load. [PR #793](https://github.com/pypsa-meets-earth/pypsa-earth/pull/793)

* Add NorthAmerican and Earth cutouts, and improve African cutout [PR #813](https://github.com/pypsa-meets-earth/pypsa-earth/pull/813)

* Bug fixing to restore Africa execution and improve performances [PR #817](https://github.com/pypsa-meets-earth/pypsa-earth/pull/817)

* Add Asian cutout [PR #826](https://github.com/pypsa-meets-earth/pypsa-earth/pull/826)

* Add a cutout for Western Asia [PR #837](https://github.com/pypsa-meets-earth/pypsa-earth/pull/837)

* Add osm_config yaml file [PR #822](https://github.com/pypsa-meets-earth/pypsa-earth/pull/822)

* Re-enable offshore wind and revise hydro [PR #830](https://github.com/pypsa-meets-earth/pypsa-earth/pull/830)

* Add databundle of cutouts for Kazakhstan for CI test  [PR #856](https://github.com/pypsa-meets-earth/pypsa-earth/pull/856). The bundle (~5MB) is used in pypsa-kz-data repository during CI tests.

* Option to specify a global upper capacity limit (using existing BAU functionality) [PR #857](https://github.com/pypsa-meets-earth/pypsa-earth/pull/857)

* Add cluster options `all`, `min` and [flex` `PR #848](https://github.com/pypsa-meets-earth/pypsa-earth/pull/857)

* Add commit id of pypsa earth in the n.meta of the .nc file per default [PR #863](https://github.com/pypsa-meets-earth/pypsa-earth/pull/863)

# PyPSA-Earth 0.2.2

**New Features and major Changes (8th July 2023)**

* Fix Natural Gas assignment bug in build_powerplants rule [PR #754](https://github.com/pypsa-meets-earth/pypsa-earth/pull/754).

* Add GEM datasets to the powerplantmatching config [PR #750](https://github.com/pypsa-meets-earth/pypsa-earth/pull/750).

* Add merge and replace functionalities when adding custom powerplants [PR #739](https://github.com/pypsa-meets-earth/pypsa-earth/pull/739). "Merge" combined the powerplantmatching data with new custom data. "Replace" allows to use fully self-collected data.

* Add functionality of attaching existing renewable caapcities from custom_powerplants.csv. [PR #744](https://github.com/pypsa-meets-earth/pypsa-earth/pull/744). If custom_powerplants are enabled and custom_powerplants.csv contains wind or solar powerplants, then p_nom and p_nom_min for renewables are extracted from custom_powerplants.csv, aggregated for each bus, and set.

* Fix dask parallel computations for e.g. cutouts calculations. Now again more than 1 core will be used when available that can lead to ~8x speed ups with 8 cores [PR #734](https://github.com/pypsa-meets-earth/pypsa-earth/pull/734) and [PR #761](https://github.com/pypsa-meets-earth/pypsa-earth/pull/761).

* Add the usage of custom rules. Custom rule files must be specified in the config as a list, e.g. custom rules: ["my_rules.smk"]. Empty by default (i.e. no custom rules). [PR #755](https://github.com/pypsa-meets-earth/pypsa-earth/pull/755)

* Add trailing whitespace linter which removes unnecessary tabs when running `pre-commit` [PR #762](https://github.com/pypsa-meets-earth/pypsa-earth/pull/762)

* Add codespell linter which corrects word spellings [PR #763](https://github.com/pypsa-meets-earth/pypsa-earth/pull/763)

* Remove RES addition functionality from attach_conventional_generators [PR #769](https://github.com/pypsa-meets-earth/pypsa-earth/pull/769). Currently wind and solar powerplants stored in powerplants.csv are added to the network by attach_conventional_generators.

* Add functionalities to download and extract emission of countries. [PR #748 https://github.com/pypsa-meets-earth/pypsa-earth/pull/748`

# PyPSA-Earth 0.2.1

**New Features and major Changes (20th May 2023)**

* Fix bug. Add graphviz to docs to compile workflows in the documentation and adapt release notes `PR #719](<https://github.com/pypsa-meets-earth/pypsa-earth/pull/719>)

* License change from GPL to AGPL as voted [here](https://github.com/pypsa-meets-earth/pypsa-earth/issues/693)

* Fix hard-coded simplification of lines to 380kV [PR #732](https://github.com/pypsa-meets-earth/pypsa-earth/pull/732).
  It is now possible to simplify the network to any other voltage level with config option [base_voltage`.

* Fix a KeyError in simplify_links caused by misinterpretation of AC lines as DC ones `PR #740](<https://github.com/pypsa-meets-earth/pypsa-earth/pull/740>).

# PyPSA-Earth 0.2.0

**New Features and major Changes (7th May 2023)**

* Finalize package restructuring [PR #462](https://github.com/pypsa-meets-earth/pypsa-earth/pull/462)

* Fix made in config.default and config.tutorial changing Monte-Carlo from true to false [PR #463](https://github.com/pypsa-meets-earth/pypsa-earth/pull/463)

* Add new config test design. It is now easy and light to test multiple configs [PR #466](https://github.com/pypsa-meets-earth/pypsa-earth/pull/466)

* Revision of documentation [PR #471](https://github.com/pypsa-meets-earth/pypsa-earth/pull/471)

* Move to new GADM version [PR #478](https://github.com/pypsa-meets-earth/pypsa-earth/pull/478)

* Update natura tiff to global scale, revise default databundle description and remove old limitations to environment [PR #470](https://github.com/pypsa-meets-earth/pypsa-earth/pull/470) and [PR #500](https://github.com/pypsa-meets-earth/pypsa-earth/pull/500)

* Update docs on installation [PR #498](https://github.com/pypsa-meets-earth/pypsa-earth/pull/498)

* Update docs on tutorial [PR #507](https://github.com/pypsa-meets-earth/pypsa-earth/pull/507)

* Moved from pycountry to country_converter [PR #493](https://github.com/pypsa-meets-earth/pypsa-earth/pull/493)

* Fix workflow in order to solve the landlock countries bug  [PR #481](https://github.com/pypsa-meets-earth/pypsa-earth/pull/481) and [PR #517](https://github.com/pypsa-meets-earth/pypsa-earth/pull/517)

* Add meta data of config to pypsa network per default. Allows keeping track of the config used to generate the network [PR #526](https://github.com/pypsa-meets-earth/pypsa-earth/pull/526)

* Fix renewable profiles generation for possible data loss in ERA5-derived cutouts [PR #511](https://github.com/pypsa-meets-earth/pypsa-earth/pull/511)

* Adapt dependencies of powerplantmatching to the PyPSA main branch [PR #527](https://github.com/pypsa-meets-earth/pypsa-earth/pull/527)

* Calculate the outputs of retrieve_databundle dynamically depending on settings [PR #529](https://github.com/pypsa-meets-earth/pypsa-earth/pull/529)

* Fix shape bug in the Voronoi cell creation [PR #541](https://github.com/pypsa-meets-earth/pypsa-earth/pull/541)

* Adapt dependencies on PyPSA to the PyPSA main branch [PR #538](https://github.com/pypsa-meets-earth/pypsa-earth/pull/538)

* Fix None geometries into regions [PR #546](https://github.com/pypsa-meets-earth/pypsa-earth/pull/546)

* Swap OpenStreetMap python download interface from esy-osm to earth-osm [PR #547](https://github.com/pypsa-meets-earth/pypsa-earth/pull/547)

* Restore saving of logger outputs [PR #559](https://github.com/pypsa-meets-earth/pypsa-earth/pull/559)

* Techno-economic parameters of technologies (e.g. costs and efficiencies) can be now retrieved from a separate repository [PyPSA/technology-data](https://github.com/pypsa/technology-data)
  that collects assumptions from a variety of sources. It is activated by default with [`enable: retrieve_cost_data: true` and controlled with `costs: year:` and `costs: version:`.
  The location of this data changed from `data/costs.csv` to `resources/costs.csv``. Adapted from [`#184](https://github.com/PyPSA/pypsa-eur/pull/184)].

* Added approaches to process contended areas [PR #572](https://github.com/pypsa-meets-earth/pypsa-earth/pull/572)

* Improve parallel capabilities of build_shapes to enable parallelization even within a country shape [PR #575](https://github.com/pypsa-meets-earth/pypsa-earth/pull/575)

* Add pypsa-eur scenario management [PR #577](https://github.com/pypsa-meets-earth/pypsa-earth/pull/577)

* Minor bug fixing and improvements [PR #580](https://github.com/pypsa-meets-earth/pypsa-earth/pull/580)

* Streamline default configuration file [PR #589](https://github.com/pypsa-meets-earth/pypsa-earth/pull/589)

* Fix rule run_test, remove code duplication, add gitstars to readme `PR #593 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/593>[

* Add new build_demand_profiles.py. It builds demand_profiles.csv and allow easier interfacing of new data `PR #582](<https://github.com/pypsa-meets-earth/pypsa-earth/pull/582>)

* Upgrade technology data to v0.5.0 [PR #600](https://github.com/pypsa-meets-earth/pypsa-earth/pull/600)

* Update simplify_network and cluster_network according to PyPSA-Eur developments [PR #597](https://github.com/pypsa-meets-earth/pypsa-earth/pull/597)

* Revise OSM cleaning to improve the cleaning process and error resilience [PR #620](https://github.com/pypsa-meets-earth/pypsa-earth/pull/620)

* Fix isolated buses when simplifying the network and add clustering by networks [PR #632](https://github.com/pypsa-meets-earth/pypsa-earth/pull/632)

* Include hydro runoff normalization [PR #631](https://github.com/pypsa-meets-earth/pypsa-earth/pull/631)

* Add REUSE compatibility [PR #651](https://github.com/pypsa-meets-earth/pypsa-earth/pull/651)

* Fix bug of missing GitHub issue template [PR #660](https://github.com/pypsa-meets-earth/pypsa-earth/pull/660)

* Fix GADM bug when using alternative clustering and store gadm shape with two letter instead of three letter ISO code  [PR #670](https://github.com/pypsa-meets-earth/pypsa-earth/pull/670)

* Fix GADM naming bug related to level-2 clustering [PR #684](https://github.com/pypsa-meets-earth/pypsa-earth/pull/684)

* Fix append bug in build_powerplants rule [PR #686](https://github.com/pypsa-meets-earth/pypsa-earth/pull/686)

* Add *zenodo_handler.py* to update and upload files via code [PR #688](https://github.com/pypsa-meets-earth/pypsa-earth/pull/688)

* Fix a few typos in docstrings [PR #695](https://github.com/pypsa-meets-earth/pypsa-earth/pull/695)

* Update and improve configuration section in documentation [PR #694](https://github.com/pypsa-meets-earth/pypsa-earth/pull/694)

* Improve earth coverage and add improve make_statistics coverage [PR #654](https://github.com/pypsa-meets-earth/pypsa-earth/pull/654)

* Fix bug for missing renewable profiles and generators [PR #714](https://github.com/pypsa-meets-earth/pypsa-earth/pull/714)

* Update instructions on how to write documentation. [PR #720](https://github.com/pypsa-meets-earth/pypsa-earth/pull/720)

* Enable workflow to run including countries with empty OSM data, test on all UN countries [PR #701 https://github.com/pypsa-meets-earth/pypsa-earth/pull/701`__

# PyPSA-Earth 0.1.0

Model rebranded from PyPSA-Africa to PyPSA-Earth. Model is part of the now called PyPSA meets Earth initiative which hosts multiple projects.

**New features and major changes (10th September 2022)**

* Identify DC lines but temporary transform them back into AC `PR #348](<https://github.com/pypsa-meets-earth/pypsa-earth/pull/348>)

* Get renewable capacities from IRENA statistics [PR #343](https://github.com/pypsa-meets-earth/pypsa-earth/pull/343)

* Bug fixing (script retrieve_databundle) and rule run_test to ease testing [PR #322](https://github.com/pypsa-meets-earth/pypsa-earth/pull/322)

* Handling non-numerical entries in raw OSM data: [PR #287](https://github.com/pypsa-meets-earth/pypsa-earth/pull/287)

* General user experience improvements: [PR #326](https://github.com/pypsa-meets-earth/pypsa-earth/pull/326)

* Fix minor validation notebook inaccuracy: [PR #332](https://github.com/pypsa-meets-earth/pypsa-earth/pull/332)

* Make clean_osm_data script work with land-locked country: [PR #341](https://github.com/pypsa-meets-earth/pypsa-earth/pull/341)

* Add demand validation notebook for 2030 prediction: [PR #344](https://github.com/pypsa-meets-earth/pypsa-earth/pull/344)

* Revise build_powerplants with new version of powerplantmatching: [PR #342](https://github.com/pypsa-meets-earth/pypsa-earth/pull/342)

* Fix typo causing the wrong coordinate reference systems (CRS) to be used when determining available land types using CLC [PR #345](https://github.com/pypsa-meets-earth/pypsa-earth/pull/345)

* Add high resolution population raster via API: [PR #325](https://github.com/pypsa-meets-earth/pypsa-earth/pull/325)

* Fix bounds of cutouts aka weather cells: [PR #347](https://github.com/pypsa-meets-earth/pypsa-earth/pull/347)

* Add new countries and update iso code: [PR #330](https://github.com/pypsa-meets-earth/pypsa-earth/pull/330)

* Fix solar pv slope and add correction factor for wake losses: [PR #335](https://github.com/pypsa-meets-earth/pypsa-earth/pull/350)

* Add renewable potential notebook: [PR #351](https://github.com/pypsa-meets-earth/pypsa-earth/pull/351)

* Make cutout workflow simpler: [PR #352](https://github.com/pypsa-meets-earth/pypsa-earth/pull/352)

* Add option to run workflow without pop and gdp raster: [PR #353](https://github.com/pypsa-meets-earth/pypsa-earth/pull/353)

* Add latitude_optimal to get optimal solar orientation by default: [Commit 1b2466b](https://github.com/pypsa-meets-earth/pypsa-earth/commit/de7d32be8807e4fc42486a60184f45680612fd46)

* Harmonize CRSs by options: [PR #356](https://github.com/pypsa-meets-earth/pypsa-earth/pull/356)

* Fix powerplantmatching problem for DRC and countries with multi-word name: [PR #359](https://github.com/pypsa-meets-earth/pypsa-earth/pull/359)

* Change default option for build_natura: [PR #360](https://github.com/pypsa-meets-earth/pypsa-earth/pull/360)

* Add renewable potential validation notebook and update others: [PR #363](https://github.com/pypsa-meets-earth/pypsa-earth/pull/363) and [PR #369](https://github.com/pypsa-meets-earth/pypsa-earth/pull/363)

* Constrain rasterio version and add plotting dependencies: [PR #365](https://github.com/pypsa-meets-earth/pypsa-earth/pull/365)

* Change solar power density form 1.7 to 4.6 MW/km2: [PR #364](https://github.com/pypsa-meets-earth/pypsa-earth/pull/364)

* Bug fixing of unexpected float value in build_powerplants: [PR #372](https://github.com/pypsa-meets-earth/pypsa-earth/pull/372) and [PR #373](https://github.com/pypsa-meets-earth/pypsa-earth/pull/373)

* Revise hydro capacities, add hydro validation notebook and minor revisions: [PR #366](https://github.com/pypsa-meets-earth/pypsa-earth/pull/366)

* Revise dropnan for regions: [PR #366](https://github.com/pypsa-meets-earth/pypsa-earth/pull/366)

* Fix bug in GADM clustering. Missing crs input: [PR #379](https://github.com/pypsa-meets-earth/pypsa-earth/pull/379)

* Optimise [availabilitymatrix` speed by factor 4-5: `PR #380](https://github.com/pypsa-meets-earth/pypsa-earth/pull/380)

* Fix bug in inline documentation for GADM and Voronoi clustering: [PR #384](https://github.com/pypsa-meets-earth/pypsa-earth/pull/384)

* Fix simple clustering enabling the creation of networks such [regions_onshore_elec_s54_14.nc`:`PR #386](https://github.com/pypsa-meets-earth/pypsa-earth/pull/386)

* Add transformer components which connect different voltage level lines: [PR #389](https://github.com/pypsa-meets-earth/pypsa-earth/pull/389)

* Enable the use of a float value for the scale in load_options: [PR #397](https://github.com/pypsa-meets-earth/pypsa-earth/pull/397)

* Add operational reserve margin according to PyPSA-Eur: [PR #399](https://github.com/pypsa-meets-earth/pypsa-earth/pull/399)

* Add optional normalization of hydro inflows by hydro_capacities or eia stats: [PR #376](https://github.com/pypsa-meets-earth/pypsa-earth/pull/376)

* Enable DC carrier in the network model and include converters into the model: [PR #392](https://github.com/pypsa-meets-earth/pypsa-earth/pull/392)

* Implement PyPSA-Eur improvements. Add gas limit constraints, add marginal cost sweeps wildcard, add and harmonize aggregation strategies, improve config usability by carrier clarifications, ease debugging by removing snakemake inputs from functions: [PR #402](https://github.com/pypsa-meets-earth/pypsa-earth/pull/402)

* Fix and add docs. Fix incomplete tutorial, recommend mamba for installation, add YouTube videos [PR #412](https://github.com/pypsa-meets-earth/pypsa-earth/pull/412) and [PR #423](https://github.com/pypsa-meets-earth/pypsa-earth/pull/423)

* Restructure the package to ease readability and fix google drive downloading method: [PR #355](https://github.com/pypsa-meets-earth/pypsa-earth/pull/355)

* Update config links to adhere to the new structure of the package: [PR #420](https://github.com/pypsa-meets-earth/pypsa-earth/pull/420)

* Improve and finalize capacity_validation notebook: [PR #406](https://github.com/pypsa-meets-earth/pypsa-earth/pull/406) and [PR #455](https://github.com/pypsa-meets-earth/pypsa-earth/pull/455)

* Fix hydro technology with the GADM clustering approach: [PR #428](https://github.com/pypsa-meets-earth/pypsa-earth/pull/428)

* Adapt for a custom shapefile for MA as a first step towards generalizing the feature: [PR #429](https://github.com/pypsa-meets-earth/pypsa-earth/pull/429)

* Improve line augmentation for network expansion explorations. Use k-edge augmenation for AC lines and random sampling for long HVDC lines: [PR #427](https://github.com/pypsa-meets-earth/pypsa-earth/pull/427)

* Fix minor bug in clustering about missing prefix assignment [PR #434](https://github.com/pypsa-meets-earth/pypsa-earth/pull/434)

* Fix major aggregation bug and adjust config: [PR #435](https://github.com/pypsa-meets-earth/pypsa-earth/pull/435)

* Fix nan techtype and wrong tech for nuclear which improves the representation of existing powerplants [PR #436](https://github.com/pypsa-meets-earth/pypsa-earth/pull/436)

* Add notebook to compare results by different solvers [PR #421](https://github.com/pypsa-meets-earth/pypsa-earth/pull/421)

* Fix overestimation of the network capacity by simplify network [PR #443](https://github.com/pypsa-meets-earth/pypsa-earth/pull/443)

* Fix output electricity column in clean_data [PR #441](https://github.com/pypsa-meets-earth/pypsa-earth/pull/441)

* Bug fixing to download global OSM and shape data: [PR #433](https://github.com/pypsa-meets-earth/pypsa-earth/pull/433)

# PyPSA-Africa 0.0.2

**New features and major changes (6th April 2022)**

* Plotting and summary features: [PR #211](https://github.com/pypsa-meets-earth/pypsa-earth/pull/211) and [PR #214](https://github.com/pypsa-meets-earth/pypsa-earth/pull/214)

* Templates for issue, PR, feature request: [PR #216](https://github.com/pypsa-meets-earth/pypsa-earth/pull/216)

* Attach hydro enabled with all hydro types: [PR #232](https://github.com/pypsa-meets-earth/pypsa-earth/pull/232)

* Parallel download of osm data: [PR #232](https://github.com/pypsa-meets-earth/pypsa-earth/pull/232)

* Decoupling iso coding from geofabrik; rule download_osm_data extended to the world: [PR #236](https://github.com/pypsa-meets-earth/pypsa-earth/pull/236)

* Rule build_shape extended to the world: [PR #236](https://github.com/pypsa-meets-earth/pypsa-earth/pull/236)

* Validation of geofabrik links: [PR #249](https://github.com/pypsa-meets-earth/pypsa-earth/pull/249)

* Generalized version of Data retrieval with google and zenodo hosting platforms: [PR #242](https://github.com/pypsa-meets-earth/pypsa-earth/pull/242) and [PR #260](https://github.com/pypsa-meets-earth/pypsa-earth/pull/260)

* Fix random state for kmean clustering, adopted from [PR 313](https://github.com/PyPSA/pypsa-eur/pull/313)

* Implement area exclusions based on land type using the Copernicus Land Cover: [PR #272](https://github.com/pypsa-meets-earth/pypsa-earth/pull/272).

* Flexible demand extraction for multiple years across the globe: [PR #275](https://github.com/pypsa-meets-earth/pypsa-earth/pull/275)

* Add CI caching and windows CI: [Commit CI windows](https://github.com/pypsa-meets-earth/pypsa-earth/commit/c98cb30e828cfda17692b8f5e1dd8e39d33766ad),  [PR #277](https://github.com/pypsa-meets-earth/pypsa-earth/pull/277).

* Change config to allow weather year extraction from snapshots as default: [PR #301](https://github.com/pypsa-meets-earth/pypsa-earth/pull/301).

* Replace Restyler by .pre-commit [PR #307 https://github.com/pypsa-meets-earth/pypsa-earth/pull/307`__.

* Solved the issue of "overpassing nodes" and restyling osm_build_network: `PR #294](<https://github.com/pypsa-meets-earth/pypsa-earth/pull/294>)

* Revise deprecations in build_shape: [PR #315](https://github.com/pypsa-meets-earth/pypsa-earth/pull/315)

# PyPSA-Africa 0.0.1

This is the first release of PyPSA-Africa which heavily builds on [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur).

**New features and major changes (24th December 2021)**

* Include new data streams for Africa model

* Demand data implementation from [GEGIS](https://github.com/pypsa-meets-earth/pypsa-earth/blob/9acf89b8756bb60d61460c1dad54625f6a67ddd5/scripts/add_electricity.py#L221-L259). Demand can be chosen for weather years and socioeconomic [ssp` scenarios

* Network is built, cleaned and processed solely on `OpenStreetMap data](<https://github.com/pypsa-meets-earth/pypsa-earth/blob/9acf89b8756bb60d61460c1dad54625f6a67ddd5/scripts/osm_pbf_power_data_extractor.py>)

* Voronoi regions, where data is aggregated towards, can be replaced by administrative [GADM zones](https://github.com/pypsa-meets-earth/pypsa-earth/commit/4aa21a29b08c4794c5e15d4209389749775a5a52)

* [Augmented line expansion feature](https://github.com/pypsa-meets-earth/pypsa-earth/pull/175) can make network meshed, connect isolated mini-grids to the main-grid.

* Community moved to [Discord](https://discord.gg/AnuJBk23FU).

* Most meeting and agenda's are [open](https://github.com/pypsa-meets-earth/pypsa-earth#get-involved).

# Release Process

* Checkout a new release branch [`git checkout -b release-v0.x.x`.

* Finalise release notes at `doc/release_notes.rst`.

* Make sure thah pinned versions of the environments `*-pinned.yaml` in `envs` folder are up-to-date.

* Update version number in `doc/conf.py`, `default.config.yaml`, `tutorial.config.yaml` and `test/config.*.yaml`.

* Open, review and merge pull request for branch `release-v0.x.x`.
  Make sure to close issues and PRs or the release milestone with it (e.g. closes #X).
  Run `pre-commit run --all`` locally and fix any issues.

* Update and checkout your local `main` and tag a release with `git tag v0.x.x`, `git push`, `git push --tags`. Include release notes in the tag message using Github UI.

* Upload code to `zenodo code repository](<https://doi.org>) with [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html).

* Create pre-built networks for [`config.default.yaml` by running `snakemake -j 1 extra_components_all_networks``.

* Upload pre-built networks to `zenodo data repository](<https://doi.org/10.5281/zenodo.3601881>) with [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) license.

* Send announcement on the [PyPSA-Earth Discord channel](https://discord.gg/AnuJBk23FU).
