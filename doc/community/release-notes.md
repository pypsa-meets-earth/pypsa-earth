


Release Notes


# Upcoming release

This part of documentation collects descriptive release notes to capture the main improvements introduced by developing the model before the next release.

**New Features and Major Changes**

* Add an option to redefine countries into subregions in [`cluster_networks`` `PR #1542](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1542)_

* Added Section covering Optimization for Energy Systems Models [PR #1558](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1558)_

* Add the configuration [`co2_budget`` to set CO₂ emission targets in multiple planning horizon years `PR #1553](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1553)_

* Add pseudo [branch()` to streamline snakemake workflow. Replace this with Snakemake implementation if Snakemake >= 8.3.0. Added an option to configure which sector components to include which when disabled, irrelevant rules are skipped automatically `PR #1538](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1538)_

* Extend functionality of isolated nodes to networks and fix to demand scaling when gdp or pop is empty [PR #1540](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1540)_

* Redesign temporal matching including yearly, monthly and hourly matching [PR #1463](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1463)_

* Add option to ignore loading network data in clean_osm_data [PR #1580](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1580)_

* Add ll wildcard option 'l{factor}' to enable line-wise transmission expansion, including optional lower and upper bounds [PR #1592](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1592)_


**Minor Changes and bug-fixing**

* Refine load shedding capacity calculation to use bus-specific maximum loads instead of fixed large values, improving solver performance and numeric stability. Rename load shedding generator carrier to 'load shedding' (now) instead of just 'load' (previously) [PR #1581](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1581)_

* The configuration setting for [`focus_weights` has been moved from `focus_weights:` to `cluster_options: focus_weights:``. Backwards compatibility to old config files is maintained `PR #1565](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1565)_

* Disable distribute_cluster only for one subnetwork and country [PR #1539](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1539)_

* Reintroduce sanitize_carriers and sanitize_location to reduce the number of warnings related to components with undefined carriers [PR #1555](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1555)_

* Fix the offwind depth calculation by providing the proper GEBCO file [PR #1559](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1559)_

* bug-fix in [`clean_osm_data.py`` so that cleaned data has the same number of circuits as in the raw data and assumptions are correctly applied `PR #1552](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1552)_

* Avoid using variable "nodes" in prepare_sector_network unless explicitly in the arguments [PR #1575](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1575)_

* Fix and reactivate the option for a custom busmap in [`cluster_network.py`` `PR #1537](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1537)_

* Bug-fixing override co2opt in add_co2_budget [PR #1597](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1597)_

* Fix shipping and aviation implementation in multi-country models [PR #1582](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1582)_

# PyPSA-Earth 0.7.0

**New Features and Major Changes**

* Revise bundle_cutouts_northamerica [PR #1479](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1479)_

* Revise pinned to locked environments [PR #1458](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1458)_

* Update to [technology-data` version v0.12.0 `PR #1452](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1452)_

* Add a command-line interface in [`databundle_cli.py` that will be triggered if `retrieve_databundle_light.py`` fails to retrieve all necessary files, providing a fallback to assist with debugging the issue `PR #1366](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1366)_

* Add Wikipedia as a source for the preparation of transport_data.csv  [PR #1410](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1410)_

* Add a new page for sector-coupled tutorial [PR #1374](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1374)_

* Add an option to redefine countries into subregions which can be activated at various stages of the workflow. The subregions can be defined either by the GADM_ID or a custom map. Currently, it is only used in the [`simplify_network`` rule `PR #1300](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1300)_

* Drop duplication of retrieve_data and COST_DIR, add params and update technology-data version [PR #1249](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1249)_

* In alternative clustering, generate hydro inflows by shape and avoid hydro inflows duplication for plants installed in the same node [PR #1120](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1120)_

* Add a function to calculate length-based efficiencies and apply it to the H2 pipelines [PR #1192](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1192)_

* Support of Linopy for Power and Sector-Coupled Modelling and latest PyPSA version [PR #1172](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1172)_

* Update workflow to geopandas >= 1.0 [PR #1276](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1276)_

* Index hydro units by their location and drop use of alternative clustering option for hydro [PR #1331](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1331)_

* Introduce universal currency conversion to allow use of currencies other than EUR in input data and results [PR #1319](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1319)_

* Introduce US-specific cost configurations and update to [technology_data_version` v0.13.2 `PR #1448](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1448)_

* Add new hydrogen production technologies and reorganize existing structure [PR #1227](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1227)_

**Minor Changes and bug-fixing**

* Remove duplicates in environment file [PR #1473](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1473)_

* Fix fallback logic for missing currency data by correctly applying the default exchange rate and build cache for storing commonly used exchange rates [PR #1492](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1492)_

* Update code to reflect the Wikipedia source page for [prepare_transport_data_input` `PR #1486](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1486)_

* Add efficiency gain and growth rates for other energy consumption and fill missing NaNs with 0 [PR #1468](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1468)_

* Fix missing nodes in prepare_sector_network [PR #1432](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1432)_

* Fix missing country attribute in cluster_network [PR #1443](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1443)_

* Fix numpy TypeErrors [PR #1430](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1430)_

* Update description of the sector-coupled capabilities [PR #1449](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1449)_

* Modified how_to_contribute documentation page [PR #1439](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1439)_

* Redefine h2export wildcard_constraint to allow digits + optional dot followed by more digites (decimals) [PR #1434](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1434)_

* Update the AL_production.csv data from another source with more countries [PR #1428](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1428)_

* Fix params for prepare_sector_network script [PR #1427](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1427)_

* Add option to adjust load shedding costs [PR #1403](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1403)_

* Fix: Use the mean value instead of the sum to remove duplicates in the urban percentage data per country and year [PR #1420](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1420)_

* Fix prepare_ports script to generate both outputs if export ports with custom data is selected [PR #1424](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1424)_

* Update the applications list for PyPSA-Earth model [PR #1413](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1413)_

* Fix problem with a discontinued World Bank data source in prepare_transport_data_input [PR #1410](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1410)_

* Fix bug in myopic run [PR #1369](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1369)_

* Fix missing focus_weights on cluster_network params without augmented_line_connections, minimize warnings on subregions [PR #1348](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1348)_

* Align structure of the components with consistency checks in the updated PyPSA version [PR #1315](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1315)_

* Prevent computation of powerplantmatching if replace option is selected for custom_powerplants [PR #1281](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1281)_

* Fix overlapping bus regions when alternative clustering is selected [PR #1287](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1287)_

* Fix readthedocs by explicitly specifying the location of the Sphinx config [PR #1292](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1292)_

* Fix lossy bidirectional links, especially H2 pipelines, which would sometimes gain H2 instead of losing it [PR #1192](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1192)_

* Fix the need for administrative rights on Windows by changing all shadow directory settings for Windows in the Snakefile [PR #1295](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1295)_ and  [PR #1301](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1301)_

* Add CITATION.cff to guide users on how cite PyPSA-Earth. Insert the DOI badge in the README linking to the very first PyPSA-Earth paper [PR #1316](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1316)_

* Remove pyomo from cluster_network and introduce scip [PR #1320](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1323)_

* Fix mixup of plant locations for hydro inflow and inflow overestimation [PR #1322](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1322)_

* Fix namibia geofk, line country tag mismatch and minor fixes [PR #1330](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1330)_

* Minor revision docker workflow to have it working on upstream [PR #1343](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1343)_

* Fix the scaling factor for time-varying loads of the sector model [PR #1372](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1372)_

* Integrate EIA data for US-specific CAGR and fuel shares for the sector model [PR #1372](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1400)_

* Revise naming of Wikipedia data for vehicles [PR #1422](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1422)_

* Monte Carlo: move qmc.discrepancy to logging only [PR #1418](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1418)_

* Extending powerplant filter option to custom powerplants [PR #1465](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1465)_

* Fix an issue with the GEBCO file by limiting libgdal-core<3.10 [PR #1519](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1519)_

* Fix a naming issue with European cutout [PR #1530](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1530)_

* Avoid adding CO2 pipeline links when option is disabled [PR #1504](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1504)_

* Add possibility to overwrite cost attributes for sector model [PR #1567](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1567)_

# PyPSA-Earth 0.6.0

**New Features and Major Changes (24th December 2024)**

* Include option in the config to allow for custom airport data [PR #1241](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1241)_

* Added Dev Containers and docker as an option to get started with pypsa-earth [PR #1228](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1228)_

* Add a list of PyPSA-Earth applications in academic and industrial projects [PR #1255](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1255)_

* Computational improvements of build_osm_network [PR #845](https://github.com/pypsa-meets-earth/pypsa-earth/pull/845)_

* Boost computational performances of set_lines_ids with cKDTree by scipy [PR #806](https://github.com/pypsa-meets-earth/pypsa-earth/pull/806)_

* Boost computational performances of set_substation_ids using DBSCAN [PR #799](https://github.com/pypsa-meets-earth/pypsa-earth/pull/799)_

* Boost computational performances of fix_overpassing_line [PR #807](https://github.com/pypsa-meets-earth/pypsa-earth/pull/807)_

**Minor Changes and bug-fixing**

* Added electricity bus to Fischer-Tropsch in prepare_sector_network.py [PR #1226](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1226)_

* Update BW, NG and BJ tutorial databundles to include gadm-like sources from geoboundaries [PR #1257](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1257)_


# PyPSA-Earth 0.5.0

**New Features and Major Changes (14th December 2024)**

* Added capabilities of cross-sectoral modeling by merging with PyPSA-Earth-Sec model [https://github.com/pypsa-meets-earth/pypsa-earth-sec`__

* The workflow configuration now supports incremental changes to the default configuration in the `config.yaml` and configfiles passed to snakemake via `--configfile myconfig.yaml`. Therefore the user may now only include settings in their `config.yaml` which differ from the default configuration. One can think of the new `config.yaml` as of a list of arguments in a python function that already have a default. So in principle the `config.yaml` could now be empty, and the workflow would still run. `PR #1053](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1053)_

* Include option of endogenous export, which optimizes the export quantity based on price signals [PR #1201](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1201)_

* Remove elec-based H2 and battery technologies before addition in [prepare_sector_network.py` script and fix bus names for links that models H2 repuspose network `PR #1198](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1198)_

* Add electricity distribution grid with solar rooftop and home battery technologies [PR #1221](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1221)_

* Include a dedicated cutout for North America in bundle_config.yaml [PR #1121](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1121)_

* Include a dedicated cutout for Europe in bundle_config.yaml [PR #1125](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1125)_

* Include a dedicated cutout for Oceania in bundle_config.yaml [PR #1157](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1157)_

* Integrate RDIR into sector rules to store intermediate data in scenario folders [PR #1154](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1154)_

* The computation of [hydro_profile.nc` in `build_renewable_profiles.py` is not differentiated whether alternative clustering is applied or not; the indexing of the different power plants in `add_electricity.py` is performed according to the bus either in case alternative clustering is applied or not and a `hydro_inflow_factor` is computed prior to the computation of `inflow_t` to split the inflow according to the capacity of each different unit of each power plant (if more units are present). `PR #1119](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1119)_

* Use BASE_DIR in rules and [_helpers.py` script for facilitate module import in subworkflow `PR #1137](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1137)_

* Enable sector rules import in subworkflow [PR #1178](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1178)_

**Minor Changes and bug-fixing**

* The default configuration for [electricity:estimate_renewable_capacities:year` was updated from 2020 to 2023. `PR #1106](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1106)_

* Fix the mismatch between buses and x, y locations while creating H2 Stores [PR #1134](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1134)_

* Enable configfile specification for mock_snakemake [PR #1135](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1135)_

* Removed duplications of devendencies in environment.yaml [PR #1128](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1128)_

* Fix pre-commit docformatter python issue. [PR #1153](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1153)_

* Drop duplicate entries in [AL_production.csv` data used in `build_industry_demand` rule `PR #1143](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1143)_

* Fix bugs in [prepare_sector_network.py` related to links with H2 buses and bug of re-addition of H2 and battery carriers in present `PR #1145](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1145)_

* Drop entries that contain non-string elements in country column of [CO2_emissions_csv` data in `prepare_transport_data_input.py` script `PR #1166](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1166)_

* Local tests are now run with [make test`. This uses a `Makefile` which runs snakemake calls with different configurations. `PR #1053](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1053)_

* Adds [Dependabot](https://docs.github.com/en/code-security/getting-started/dependabot-quickstart-guide)_ to keep GitHub actions up to date. [PR #1184](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1184)_

* Adds code security scans via [CodeQL](https://codeql.github.com/)_ to CI. [PR #1185](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1185)_

* Adds CI to update keep pinned environment files up to date. [PR #1183](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1183)_ and  [PR #1210](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1210)_

* Revise ports data for export in [add_export.py` related to sector model `PR #1175](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1175)_

* Restore string values of tech_colors in config file [PR #1205](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1205)_

* Drop vrestil dependency [PR #1220](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1220)_

* Include a configuration option to simplify / not simplify shapefiles based on a boolean value specified under [build_shape_options:simplify_gadm` option in the config file `PR 1138](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1138)

* Fix the mismatch between buses and x, y locations while creating H2 Stores [PR #1134](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1134)

* Remove duplicate entries from hydrogen export ports [PR #1233](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1233)_

* Fix the environment placing a version limit to numpoly [PR #1237](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1237)_

# PyPSA-Earth 0.4.1

**New Features and Major Changes (19th September 2024)**

* Add functionality to modify the cost assumptions using config files [PR #1097](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1097)_

**Minor Changes and bug-fixing**

* Remove unused [countries_codes` argument from `load_GDP` function in `build_shapes.py` script, which was not being called as intended with positional arguments `PR #1069](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1069)_
* Fixed problematic float parsing ([_parse_float`) in `clean_osm_data.py` to make sure all OSM lines are correctly accounted for `PR #1089](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1089)_
* Fix minor bug for advanced csp implementation [PR #1076](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1076)_
* Fix minor bug in [build_powerplants.py` where the gas technology assignment incorrectly introduced NaN values for all powerplant technologies. `PR #1102](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1102)_


# PyPSA-Earth 0.4.0

**New Features and Major Changes (27th July 2024)**

* Improve Monte Carlo feature with more distributions types, independent by PyPSA component. [PR #930](https://github.com/pypsa-meets-earth/pypsa-earth/pull/930)_

* Introduce flexible regional selection of the demand files of GEGIS. [PR #991](https://github.com/pypsa-meets-earth/pypsa-earth/pull/991)_

* Generalize line types for AC and DC networks. [PR #999](https://github.com/pypsa-meets-earth/pypsa-earth/pull/999)_

* Add an option to merge isolated networks into respective backbone networks by countries. [PR #903](https://github.com/pypsa-meets-earth/pypsa-earth/pull/903)_

* Add an option to use csv format for custom demand imports. [PR #995](https://github.com/pypsa-meets-earth/pypsa-earth/pull/995)_


**Minor Changes and bug-fixing**

* Minor bug-fixing to run the cluster wildcard min [PR #1019](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1019)_

* Add option to adjust load scale for each individual countries [PR #1006](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1006)_

* Minor bug-fixing to get the generalised line types work for DC lines and AC lines. [PR #1008](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1008)_ , [PR #1011](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1011)_ and [PR #1013](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1013)_

* Minor bug-fixing for GADM_ID format naming. [PR #980](https://github.com/pypsa-meets-earth/pypsa-earth/pull/980)_, [PR #986](https://github.com/pypsa-meets-earth/pypsa-earth/pull/986)_ and [PR #989](https://github.com/pypsa-meets-earth/pypsa-earth/pull/989)_

* Fix download_osm_data compatibility for earth-osm v2.1. [PR #954](https://github.com/pypsa-meets-earth/pypsa-earth/pull/954)_ and [PR #988](https://github.com/pypsa-meets-earth/pypsa-earth/pull/988)_

* Improve geometry filtering in clean_osm_data. [PR #989](https://github.com/pypsa-meets-earth/pypsa-earth/pull/989)_

* Revise bus region definition by gadm. [PR #1001](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1001)_

* Documentation improvements. [PR #1007](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1007)_

* Remove unnecessary imports. [PR #1020](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1020)_

* Resolve pandas deprecation warning. [PR #1023](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1023)_

* Create files where the code outputs the value of the objective function. [PR #1033](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1033)_

* Introduce versioning of the configuration files. [PR #1058](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1058)_

* Fix bug for hydro inflow normalization for gadm regions (alternative clustering). [PR #1057](https://github.com/pypsa-meets-earth/pypsa-earth/pull/1057)_

* Minor bug-fixing for s_nom_min. [PR #961](https://github.com/pypsa-meets-earth/pypsa-earth/pull/961)_


# PyPSA-Earth 0.3.0

**New Features and major Changes (24th December 2023)**

* Keep all traceback in logs. [PR #898](https://github.com/pypsa-meets-earth/pypsa-earth/pull/898)_

* Function added in clean_osm_data script to allow the use of custom network data instead or on-top of OSM data. [PR #842]('https://github.com/pypsa-meets-earth/pypsa-earth/pull/842)_

* Improve retrieve_databundle to prioritize smallest databundles [PR #911](https://github.com/pypsa-meets-earth/pypsa-earth/pull/911)_

* Add functionality to load shapefiles for hydrobasins directly from the data source directly [PR #919](https://github.com/pypsa-meets-earth/pypsa-earth/pull/919)_

* Use [new CC0 v1 dataset](https://doi.org/10.7910/DVN/XIV9BL)_ for the natura input and automate download of WDPA protected planet data [PR #913](https://github.com/pypsa-meets-earth/pypsa-earth/pull/913)_

**Minor Changes and bug-fixing**

* Revise databundles and improve logging in retrieve_databundle [PR #928](https://github.com/pypsa-meets-earth/pypsa-earth/pull/928)_

* Improve documentation on installation and short tutorial [PR #918](https://github.com/pypsa-meets-earth/pypsa-earth/pull/918)_

# PyPSA-Earth 0.2.3

**New Features and major Changes (19th October 2023)**

* Add params: section in rule definition to keep track of changed settings in config.yaml. [PR #823](https://github.com/pypsa-meets-earth/pypsa-earth/pull/823)_ and [PR #880](https://github.com/pypsa-meets-earth/pypsa-earth/pull/880)_

* Fix Natural Gas implementation in "add_electricity" to avoid "Natural Gas" to be filtered out [PR #797](https://github.com/pypsa-meets-earth/pypsa-earth/pull/797)_

* Improve network simplification routine to account for representation HVDC as Line component [PR #743](https://github.com/pypsa-meets-earth/pypsa-earth/pull/743)_

* Remove deprecated pypsa.networkclustering approach and replace by pypsa.clustering.spatial [PR #786](https://github.com/pypsa-meets-earth/pypsa-earth/pull/786)_

* Drop code-dependency from vresutil [PR #803](https://github.com/pypsa-meets-earth/pypsa-earth/pull/803)_

* Add a check to ensure match between a cutout and a modelled area [PR #805](https://github.com/pypsa-meets-earth/pypsa-earth/pull/805)_

* Support renewables or renewable expansion to meet a desired share of total load. [PR #793](https://github.com/pypsa-meets-earth/pypsa-earth/pull/793)_

* Add NorthAmerican and Earth cutouts, and improve African cutout [PR #813](https://github.com/pypsa-meets-earth/pypsa-earth/pull/813)_

* Bug fixing to restore Africa execution and improve performances [PR #817](https://github.com/pypsa-meets-earth/pypsa-earth/pull/817)_

* Add Asian cutout [PR #826](https://github.com/pypsa-meets-earth/pypsa-earth/pull/826)_

* Add a cutout for Western Asia [PR #837](https://github.com/pypsa-meets-earth/pypsa-earth/pull/837)_

* Add osm_config yaml file [PR #822](https://github.com/pypsa-meets-earth/pypsa-earth/pull/822)_

* Re-enable offshore wind and revise hydro [PR #830](https://github.com/pypsa-meets-earth/pypsa-earth/pull/830)_

* Add databundle of cutouts for Kazakhstan for CI test  [PR #856](https://github.com/pypsa-meets-earth/pypsa-earth/pull/856)_. The bundle (~5MB) is used in pypsa-kz-data repository during CI tests.

* Option to specify a global upper capacity limit (using existing BAU functionality) [PR #857](https://github.com/pypsa-meets-earth/pypsa-earth/pull/857)_

* Add cluster options [all`, `min` and `flex` `PR #848](https://github.com/pypsa-meets-earth/pypsa-earth/pull/857)_

* Add commit id of pypsa earth in the n.meta of the .nc file per default [PR #863](https://github.com/pypsa-meets-earth/pypsa-earth/pull/863)_

# PyPSA-Earth 0.2.2

**New Features and major Changes (8th July 2023)**

* Fix Natural Gas assignment bug in build_powerplants rule [PR #754](https://github.com/pypsa-meets-earth/pypsa-earth/pull/754)_.

* Add GEM datasets to the powerplantmatching config [PR #750](https://github.com/pypsa-meets-earth/pypsa-earth/pull/750)_.

* Add merge and replace functionalities when adding custom powerplants [PR #739](https://github.com/pypsa-meets-earth/pypsa-earth/pull/739)_. "Merge" combined the powerplantmatching data with new custom data. "Replace" allows to use fully self-collected data.

* Add functionality of attaching existing renewable caapcities from custom_powerplants.csv. [PR #744](https://github.com/pypsa-meets-earth/pypsa-earth/pull/744)_. If custom_powerplants are enabled and custom_powerplants.csv contains wind or solar powerplants, then p_nom and p_nom_min for renewables are extracted from custom_powerplants.csv, aggregated for each bus, and set.

* Fix dask parallel computations for e.g. cutouts calculations. Now again more than 1 core will be used when available that can lead to ~8x speed ups with 8 cores [PR #734](https://github.com/pypsa-meets-earth/pypsa-earth/pull/734)_ and [PR #761](https://github.com/pypsa-meets-earth/pypsa-earth/pull/761)_.

* Add the usage of custom rules. Custom rule files must be specified in the config as a list, e.g. custom rules: ["my_rules.smk"]. Empty by default (i.e. no custom rules). [PR #755](https://github.com/pypsa-meets-earth/pypsa-earth/pull/755)_

* Add trailing whitespace linter which removes unnecessary tabs when running [pre-commit` `PR #762](https://github.com/pypsa-meets-earth/pypsa-earth/pull/762)_

* Add codespell linter which corrects word spellings [PR #763](https://github.com/pypsa-meets-earth/pypsa-earth/pull/763)_

* Remove RES addition functionality from attach_conventional_generators [PR #769](https://github.com/pypsa-meets-earth/pypsa-earth/pull/769)_. Currently wind and solar powerplants stored in powerplants.csv are added to the network by attach_conventional_generators.

* Add functionalities to download and extract emission of countries. [PR #748 https://github.com/pypsa-meets-earth/pypsa-earth/pull/748`

# PyPSA-Earth 0.2.1

**New Features and major Changes (20th May 2023)**

* Fix bug. Add graphviz to docs to compile workflows in the documentation and adapt release notes `PR #719](https://github.com/pypsa-meets-earth/pypsa-earth/pull/719)_

* License change from GPL to AGPL as voted [here](https://github.com/pypsa-meets-earth/pypsa-earth/issues/693)_

* Fix hard-coded simplification of lines to 380kV [PR #732](https://github.com/pypsa-meets-earth/pypsa-earth/pull/732)_.
  It is now possible to simplify the network to any other voltage level with config option [base_voltage`.

* Fix a KeyError in simplify_links caused by misinterpretation of AC lines as DC ones `PR #740](https://github.com/pypsa-meets-earth/pypsa-earth/pull/740)_.

# PyPSA-Earth 0.2.0

**New Features and major Changes (7th May 2023)**

* Finalize package restructuring [PR #462](https://github.com/pypsa-meets-earth/pypsa-earth/pull/462)_

* Fix made in config.default and config.tutorial changing Monte-Carlo from true to false [PR #463](https://github.com/pypsa-meets-earth/pypsa-earth/pull/463)_

* Add new config test design. It is now easy and light to test multiple configs [PR #466](https://github.com/pypsa-meets-earth/pypsa-earth/pull/466)_

* Revision of documentation [PR #471](https://github.com/pypsa-meets-earth/pypsa-earth/pull/471)_

* Move to new GADM version [PR #478](https://github.com/pypsa-meets-earth/pypsa-earth/pull/478)_

* Update natura tiff to global scale, revise default databundle description and remove old limitations to environment [PR #470](https://github.com/pypsa-meets-earth/pypsa-earth/pull/470)_ and [PR #500](https://github.com/pypsa-meets-earth/pypsa-earth/pull/500)_

* Update docs on installation [PR #498](https://github.com/pypsa-meets-earth/pypsa-earth/pull/498)_

* Update docs on tutorial [PR #507](https://github.com/pypsa-meets-earth/pypsa-earth/pull/507)_

* Moved from pycountry to country_converter [PR #493](https://github.com/pypsa-meets-earth/pypsa-earth/pull/493)_

* Fix workflow in order to solve the landlock countries bug  [PR #481](https://github.com/pypsa-meets-earth/pypsa-earth/pull/481)_ and [PR #517](https://github.com/pypsa-meets-earth/pypsa-earth/pull/517)_

* Add meta data of config to pypsa network per default. Allows keeping track of the config used to generate the network [PR #526](https://github.com/pypsa-meets-earth/pypsa-earth/pull/526)_

* Fix renewable profiles generation for possible data loss in ERA5-derived cutouts [PR #511](https://github.com/pypsa-meets-earth/pypsa-earth/pull/511)_

* Adapt dependencies of powerplantmatching to the PyPSA main branch [PR #527](https://github.com/pypsa-meets-earth/pypsa-earth/pull/527)_

* Calculate the outputs of retrieve_databundle dynamically depending on settings [PR #529](https://github.com/pypsa-meets-earth/pypsa-earth/pull/529)_

* Fix shape bug in the Voronoi cell creation [PR #541](https://github.com/pypsa-meets-earth/pypsa-earth/pull/541)_

* Adapt dependencies on PyPSA to the PyPSA main branch [PR #538](https://github.com/pypsa-meets-earth/pypsa-earth/pull/538)_

* Fix None geometries into regions [PR #546](https://github.com/pypsa-meets-earth/pypsa-earth/pull/546)_

* Swap OpenStreetMap python download interface from esy-osm to earth-osm [PR #547](https://github.com/pypsa-meets-earth/pypsa-earth/pull/547)_

* Restore saving of logger outputs [PR #559](https://github.com/pypsa-meets-earth/pypsa-earth/pull/559)_

* Techno-economic parameters of technologies (e.g. costs and efficiencies) can be now retrieved from a separate repository [PyPSA/technology-data](https://github.com/pypsa/technology-data)
  that collects assumptions from a variety of sources. It is activated by default with [`enable: retrieve_cost_data: true` and controlled with `costs: year:` and `costs: version:`.
  The location of this data changed from `data/costs.csv` to `resources/costs.csv``. Adapted from [`#184](https://github.com/PyPSA/pypsa-eur/pull/184)].

* Added approaches to process contended areas [PR #572](https://github.com/pypsa-meets-earth/pypsa-earth/pull/572)_

* Improve parallel capabilities of build_shapes to enable parallelization even within a country shape [PR #575](https://github.com/pypsa-meets-earth/pypsa-earth/pull/575)_

* Add pypsa-eur scenario management [PR #577](https://github.com/pypsa-meets-earth/pypsa-earth/pull/577)_

* Minor bug fixing and improvements [PR #580](https://github.com/pypsa-meets-earth/pypsa-earth/pull/580)_

* Streamline default configuration file [PR #589](https://github.com/pypsa-meets-earth/pypsa-earth/pull/589)_

* Fix rule run_test, remove code duplication, add gitstars to readme `PR #593 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/593>[

* Add new build_demand_profiles.py. It builds demand_profiles.csv and allow easier interfacing of new data `PR #582](https://github.com/pypsa-meets-earth/pypsa-earth/pull/582)_

* Upgrade technology data to v0.5.0 [PR #600](https://github.com/pypsa-meets-earth/pypsa-earth/pull/600)_

* Update simplify_network and cluster_network according to PyPSA-Eur developments [PR #597](https://github.com/pypsa-meets-earth/pypsa-earth/pull/597)_

* Revise OSM cleaning to improve the cleaning process and error resilience [PR #620](https://github.com/pypsa-meets-earth/pypsa-earth/pull/620)_

* Fix isolated buses when simplifying the network and add clustering by networks [PR #632](https://github.com/pypsa-meets-earth/pypsa-earth/pull/632)_

* Include hydro runoff normalization [PR #631](https://github.com/pypsa-meets-earth/pypsa-earth/pull/631)_

* Add REUSE compatibility [PR #651](https://github.com/pypsa-meets-earth/pypsa-earth/pull/651)_

* Fix bug of missing GitHub issue template [PR #660](https://github.com/pypsa-meets-earth/pypsa-earth/pull/660)_

* Fix GADM bug when using alternative clustering and store gadm shape with two letter instead of three letter ISO code  [PR #670](https://github.com/pypsa-meets-earth/pypsa-earth/pull/670)_

* Fix GADM naming bug related to level-2 clustering [PR #684](https://github.com/pypsa-meets-earth/pypsa-earth/pull/684)_

* Fix append bug in build_powerplants rule [PR #686](https://github.com/pypsa-meets-earth/pypsa-earth/pull/686)_

* Add *zenodo_handler.py* to update and upload files via code [PR #688](https://github.com/pypsa-meets-earth/pypsa-earth/pull/688)_

* Fix a few typos in docstrings [PR #695](https://github.com/pypsa-meets-earth/pypsa-earth/pull/695)_

* Update and improve configuration section in documentation [PR #694](https://github.com/pypsa-meets-earth/pypsa-earth/pull/694)_

* Improve earth coverage and add improve make_statistics coverage [PR #654](https://github.com/pypsa-meets-earth/pypsa-earth/pull/654)_

* Fix bug for missing renewable profiles and generators [PR #714](https://github.com/pypsa-meets-earth/pypsa-earth/pull/714)_

* Update instructions on how to write documentation. [PR #720](https://github.com/pypsa-meets-earth/pypsa-earth/pull/720)_

* Enable workflow to run including countries with empty OSM data, test on all UN countries [PR #701 https://github.com/pypsa-meets-earth/pypsa-earth/pull/701`__

# PyPSA-Earth 0.1.0

Model rebranded from PyPSA-Africa to PyPSA-Earth. Model is part of the now called PyPSA meets Earth initiative which hosts multiple projects.

**New features and major changes (10th September 2022)**

* Identify DC lines but temporary transform them back into AC `PR #348](https://github.com/pypsa-meets-earth/pypsa-earth/pull/348)_

* Get renewable capacities from IRENA statistics [PR #343](https://github.com/pypsa-meets-earth/pypsa-earth/pull/343)_

* Bug fixing (script retrieve_databundle) and rule run_test to ease testing [PR #322](https://github.com/pypsa-meets-earth/pypsa-earth/pull/322)_

* Handling non-numerical entries in raw OSM data: [PR #287](https://github.com/pypsa-meets-earth/pypsa-earth/pull/287)_

* General user experience improvements: [PR #326](https://github.com/pypsa-meets-earth/pypsa-earth/pull/326)_

* Fix minor validation notebook inaccuracy: [PR #332](https://github.com/pypsa-meets-earth/pypsa-earth/pull/332)_

* Make clean_osm_data script work with land-locked country: [PR #341](https://github.com/pypsa-meets-earth/pypsa-earth/pull/341)_

* Add demand validation notebook for 2030 prediction: [PR #344](https://github.com/pypsa-meets-earth/pypsa-earth/pull/344)_

* Revise build_powerplants with new version of powerplantmatching: [PR #342](https://github.com/pypsa-meets-earth/pypsa-earth/pull/342)_

* Fix typo causing the wrong coordinate reference systems (CRS) to be used when determining available land types using CLC [PR #345](https://github.com/pypsa-meets-earth/pypsa-earth/pull/345)_

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

* Fix simple clustering enabling the creation of networks such [regions_onshore_elec_s54_14.nc`: `PR #386](https://github.com/pypsa-meets-earth/pypsa-earth/pull/386)

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

* Plotting and summary features: [PR #211](https://github.com/pypsa-meets-earth/pypsa-earth/pull/211)_ and [PR #214](https://github.com/pypsa-meets-earth/pypsa-earth/pull/214)_

* Templates for issue, PR, feature request: [PR #216](https://github.com/pypsa-meets-earth/pypsa-earth/pull/216)_

* Attach hydro enabled with all hydro types: [PR #232](https://github.com/pypsa-meets-earth/pypsa-earth/pull/232)_

* Parallel download of osm data: [PR #232](https://github.com/pypsa-meets-earth/pypsa-earth/pull/232)_

* Decoupling iso coding from geofabrik; rule download_osm_data extended to the world: [PR #236](https://github.com/pypsa-meets-earth/pypsa-earth/pull/236)_

* Rule build_shape extended to the world: [PR #236](https://github.com/pypsa-meets-earth/pypsa-earth/pull/236)_

* Validation of geofabrik links: [PR #249](https://github.com/pypsa-meets-earth/pypsa-earth/pull/249)_

* Generalized version of Data retrieval with google and zenodo hosting platforms: [PR #242](https://github.com/pypsa-meets-earth/pypsa-earth/pull/242)_ and [PR #260](https://github.com/pypsa-meets-earth/pypsa-earth/pull/260)_

* Fix random state for kmean clustering, adopted from [PR 313](https://github.com/PyPSA/pypsa-eur/pull/313)_

* Implement area exclusions based on land type using the Copernicus Land Cover: [PR #272](https://github.com/pypsa-meets-earth/pypsa-earth/pull/272)_.

* Flexible demand extraction for multiple years across the globe: [PR #275](https://github.com/pypsa-meets-earth/pypsa-earth/pull/275)

* Add CI caching and windows CI: [Commit CI windows](https://github.com/pypsa-meets-earth/pypsa-earth/commit/c98cb30e828cfda17692b8f5e1dd8e39d33766ad)_,  [PR #277](https://github.com/pypsa-meets-earth/pypsa-earth/pull/277)_.

* Change config to allow weather year extraction from snapshots as default: [PR #301](https://github.com/pypsa-meets-earth/pypsa-earth/pull/301)_.

* Replace Restyler by .pre-commit [PR #307 https://github.com/pypsa-meets-earth/pypsa-earth/pull/307`__.

* Solved the issue of "overpassing nodes" and restyling osm_build_network: `PR #294](https://github.com/pypsa-meets-earth/pypsa-earth/pull/294)_

* Revise deprecations in build_shape: [PR #315](https://github.com/pypsa-meets-earth/pypsa-earth/pull/315)_


# PyPSA-Africa 0.0.1

This is the first release of PyPSA-Africa which heavily builds on [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur)_.

**New features and major changes (24th December 2021)**

* Include new data streams for Africa model

* Demand data implementation from [GEGIS](https://github.com/pypsa-meets-earth/pypsa-earth/blob/9acf89b8756bb60d61460c1dad54625f6a67ddd5/scripts/add_electricity.py#L221-L259)_. Demand can be chosen for weather years and socioeconomic [ssp` scenarios

* Network is built, cleaned and processed solely on `OpenStreetMap data](https://github.com/pypsa-meets-earth/pypsa-earth/blob/9acf89b8756bb60d61460c1dad54625f6a67ddd5/scripts/osm_pbf_power_data_extractor.py)_

* Voronoi regions, where data is aggregated towards, can be replaced by administrative [GADM zones](https://github.com/pypsa-meets-earth/pypsa-earth/commit/4aa21a29b08c4794c5e15d4209389749775a5a52)_

* [Augmented line expansion feature](https://github.com/pypsa-meets-earth/pypsa-earth/pull/175)_ can make network meshed, connect isolated mini-grids to the main-grid.

* Community moved to [Discord](https://discord.gg/AnuJBk23FU)_.

* Most meeting and agenda's are [open](https://github.com/pypsa-meets-earth/pypsa-earth#get-involved)_.


# Release Process

* Checkout a new release branch [`git checkout -b release-v0.x.x`.

* Finalise release notes at `doc/release_notes.rst`.

* Make sure thah pinned versions of the environments `*-pinned.yaml` in `envs` folder are up-to-date.

* Update version number in `doc/conf.py`, `default.config.yaml`, `tutorial.config.yaml` and `test/config.*.yaml`.

* Open, review and merge pull request for branch `release-v0.x.x`.
  Make sure to close issues and PRs or the release milestone with it (e.g. closes #X).
  Run `pre-commit run --all`` locally and fix any issues.

* Update and checkout your local `main` and tag a release with `git tag v0.x.x`, `git push`, `git push --tags`. Include release notes in the tag message using Github UI.

* Upload code to `zenodo code repository](https://doi.org) with [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html).

* Create pre-built networks for [`config.default.yaml` by running `snakemake -j 1 extra_components_all_networks``.

* Upload pre-built networks to `zenodo data repository](https://doi.org/10.5281/zenodo.3601881) with [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) license.

* Send announcement on the [PyPSA-Earth Discord channel](https://discord.gg/AnuJBk23FU).

