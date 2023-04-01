.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

##########################################
Release Notes
##########################################


Upcoming Release
================

**New Features and major Changes**

* Finalize package restructuring `PR #462 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/462>`__

* Fix made in config.default and config.tutorial changing Monte-Carlo from true to false `PR #463 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/463>`__

* Add new config test design. It is now easy and light to test multiple configs `PR #466 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/466>`__

* Revision of documenation `PR #471 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/471>`__

* Move to new GADM version `PR #478 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/478>`__

* Update natura tiff to global scale, revise default databundle description and remove old limitations to environment `PR #470 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/470>`__ and `PR #500 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/500>`__

* Update docs on installation `PR #498 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/498>`__

* Update docs on tutorial `PR #507 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/507>`__

* Moved from pycountry to country_converter `PR #493 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/493>`__

* Fix workflow in order to solve the landlock countries bug  `PR #481 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/481>`__ and `PR #517 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/517>`__

* Add meta data of config to pypsa network per default. Allows keeping track of the config used to generate the network `PR #526 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/526>`__

* Fix renewable profiles generation for possible data loss in ERA5-derived cutouts `PR #511 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/511>`__

* Adapt dependencies of powerplantmatching to the PyPSA main branch `PR #527 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/527>`__

* Calculate the outputs of retrieve_databundle dynamically depending on settings `PR #529 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/529>`__

* Fix shape bug in the Voronoi cell creation `PR #541 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/541>`__

* Adapt dependencies on PyPSA to the PyPSA main branch `PR #538 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/538>`__

* Fix None geometries into regions `PR #546 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/546>`__

* Swap OpenStreetMap python download interface from esy-osm to earth-osm `PR #547 https://github.com/pypsa-meets-earth/pypsa-earth/pull/547>`__

* Restore saving of logger outputs `PR #559 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/559>`__

* Techno-economic parameters of technologies (e.g. costs and efficiencies) can be now retrieved from a separate repository `PyPSA/technology-data <https://github.com/pypsa/technology-data>`_ 
  that collects assumptions from a variety of sources. It is activated by default with ``enable: retrieve_cost_data: true`` and controlled with ``costs: year:`` and ``costs: version:``. 
  The location of this data changed from ``data/costs.csv`` to ``resources/costs.csv``. Adapted from [`#184 <https://github.com/PyPSA/pypsa-eur/pull/184>`_].

* Added approaches to process contended areas `PR #572 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/572>`__

* Improve parallel capabilities of build_shapes to enable parallelization even within a country shape `PR #575 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/575>`__

* Add pypsa-eur scenario management `PR #577 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/577>`__

* Minor bug fixing and improvements `PR #580 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/580>`__

* Streamline default configuration file `PR #589 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/589>`__

* Fix rule run_test, remove code duplication, add gitstars to readme `PR #593 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/593>`

* Add new build_demand_profiles.py. It builds demand_profiles.csv and allow easier interfacing of new data `PR #582 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/582>`__

* Upgrade technology data to v0.5.0 `PR #600 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/600>`__

* Update simplify_network and cluster_network according to PyPSA-Eur developments `PR #597 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/597>`__

* Revise OSM cleaning to improve the cleaning process and error resilience `PR #620 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/620>`__

* Fix isolated buses when simplifying the network and add clustering by networks `PR #632 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/632>`__

* Include hydro runoff normalization `PR #631 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/631>`__

* Add REUSE compatability `PR #651 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/651>`__

* Fix bug of missing GitHub issue template `PR #660 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/660>`__

PyPSA-Earth 0.1.0
=================

Model rebranded from PyPSA-Africa to PyPSA-Earth. Model is part of the now called PyPSA meets Earth initiative which hosts mutliple projects.

**New features and major changes (10th September 2022)**

* Identify DC lines but temporary transform them back into AC `PR #348 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/348>`__

* Get renewable capacities from IRENA statistics `PR #343 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/343>`__

* Bug fixing (script retrieve_databundle) and rule run_test to ease testing `PR #322 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/322>`__

* Handling non-numerical entries in raw OSM data: `PR #287 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/287>`__

* General user experience improvements: `PR #326 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/326>`__

* Fix minor validation notebook inaccuracy: `PR #332 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/332>`__

* Make clean_osm_data script work with land-locked country: `PR #341 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/341>`__

* Add demand validation notebook for 2030 prediction: `PR #344 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/344>`__

* Revise build_powerplants with new version of powerplantmatching: `PR #342 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/342>`__

* Fix typo causing the wrong coordinate reference systems (CRS) to be used when determining available land types using CLC `PR #345 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/345>`__

* Add high resolution population raster via API: `PR #325 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/325>`_

* Fix bounds of cutouts aka weather cells: `PR #347 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/347>`_

* Add new countries and update iso code: `PR #330 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/330>`_

* Fix solar pv slope and add correction factor for wake losses: `PR #335 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/350>`_

* Add renewable potential notebook: `PR #351 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/351>`_

* Make cutout workflow simpler: `PR #352 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/352>`_

* Add option to run workflow without pop and gdp raster: `PR #353 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/353>`_

* Add latitude_optimal to get optimal solar orientation by default: `Commit 1b2466b <https://github.com/pypsa-meets-earth/pypsa-earth/commit/de7d32be8807e4fc42486a60184f45680612fd46>`_

* Harmonize CRSs by options: `PR #356 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/356>`_

* Fix powerplantmatching problem for DRC and countries with multi-word name: `PR #359 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/359>`_

* Change default option for build_natura: `PR #360 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/360>`_

* Add renewable potential validation notebook and update others: `PR #363 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/363>`_ and `PR #369 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/363>`_

* Constrain rasterio version and add plotting dependencies: `PR #365 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/365>`_

* Change solar power density form 1.7 to 4.6 MW/km2: `PR #364 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/364>`_

* Bug fixing of unexpected float value in build_powerplants: `PR #372 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/372>`_ and `PR #373 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/373>`_

* Revise hydro capacities, add hydro validation notebook and minor revisions: `PR #366 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/366>`_

* Revise dropnan for regions: `PR #366 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/366>`_

* Fix bug in GADM clustering. Missing crs input: `PR #379 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/379>`_

* Optimise `availabilitymatrix` speed by factor 4-5: `PR #380 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/380>`_

* Fix bug in inline documentation for GADM and Voronoi clustering: `PR #384 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/384>`_

* Fix simple clustering enabling the creation of networks such `regions_onshore_elec_s54_14.nc`: `PR #386 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/386>`_

* Add transformer components which connect different voltage level lines: `PR #389 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/389>`_

* Enable the use of a float value for the scale in load_options: `PR #397 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/397>`_

* Add operational reserve margin according to PyPSA-Eur: `PR #399 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/399>`_

* Add optional normalization of hydro inflows by hydro_capacities or eia stats: `PR #376 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/376>`_

* Enable DC carrier in the network model and include converters into the model: `PR #392 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/392>`_

* Implement PyPSA-Eur improvements. Add gas limit constraints, add marginal cost sweeps wildcard, add and harmonize aggregation strategies, improve config usability by carrier clarifications, ease debugging by removing snakemake inputs from functions: `PR #402 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/402>`_

* Fix and add docs. Fix incomplete tutorial, recommend mamba for installation, add YouTube videos `PR #412 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/412>`_ and `PR #423 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/423>`_

* Restructure the package to ease readability and fix google drive downloading method: `PR #355 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/355>`_

* Update config links to adhere to the new structure of the package: `PR #420 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/420>`_

* Improve and finalize capacity_validation notebook: `PR #406 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/406>`_ and `PR #455 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/455>`_

* Fix hydro technology with the GADM clustering approach: `PR #428 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/428>`_

* Adapt for a custom shapefile for MA as a first step towards generalizing the feature: `PR #429 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/429>`_

* Improve line augmentation for network expansion explorations. Use k-edge augmenation for AC lines and random sampling for long HVDC lines: `PR #427 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/427>`_

* Fix minor bug in clustering about missing prefix assignment `PR #434 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/434>`_

* Fix major aggregation bug and adjust config: `PR #435 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/435>`_

* Fix nan techtype and wrong tech for nuclear which improves the representation of existing powerplants `PR #436 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/436>`_

* Add notebook to compare results by different solvers `PR #421 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/421>`_

* Fix overestimation of the network capacity by simplify network `PR #443 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/443>`_

* Fix output electricity column in clean_data `PR #441 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/441>`_

* Bug fixing to download global OSM and shape data: `PR #433 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/433>`_

PyPSA-Africa 0.0.2 
==================

**New features and major changes (6th April 2022)**

* Plotting and summary features: `PR #211 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/211>`__ and `PR #214 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/214>`__

* Templates for issue, PR, feature request: `PR #216 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/216>`__

* Attach hydro enabled with all hydro types: `PR #232 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/232>`__

* Parallel download of osm data: `PR #232 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/232>`__

* Decoupling iso coding from geofabrik; rule download_osm_data extended to the world: `PR #236 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/236>`__

* Rule build_shape extended to the world: `PR #236 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/236>`__

* Validation of geofabrik links: `PR #249 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/249>`__

* Generalized version of Data retrieval with google and zenodo hosting platforms: `PR #242 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/242>`__ and `PR #260 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/260>`__

* Fix random state for kmean clustering, adopted from `PR 313 <https://github.com/PyPSA/pypsa-eur/pull/313>`__

* Implement area exclusions based on land type using the Copernicus Land Cover: `PR #272 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/272>`__.

* Flexible demand extraction for multiple years across the globe: `PR #275 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/275>`_

* Add CI caching and windows CI: `Commit CI windows <https://github.com/pypsa-meets-earth/pypsa-earth/commit/c98cb30e828cfda17692b8f5e1dd8e39d33766ad>`__,  `PR #277 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/277>`__.

* Change config to allow weather year extraction from snapshots as default: `PR #301 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/301>`__.

* Replace Restyler by .pre-commit `PR #307 https://github.com/pypsa-meets-earth/pypsa-earth/pull/307`__.

* Solved the issue of "overpassing nodes" and restyling osm_build_network: `PR #294 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/294>`__

* Revise deprecations in build_shape: `PR #315 <https://github.com/pypsa-meets-earth/pypsa-earth/pull/315>`__


PyPSA-Africa 0.0.1 
==================

This is the first release of PyPSA-Africa which heavily builds on `PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`__.

**New features and major changes (24th December 2021)**

* Include new data streams for Africa model

* Demand data implementation from `GEGIS <https://github.com/pypsa-meets-earth/pypsa-earth/blob/9acf89b8756bb60d61460c1dad54625f6a67ddd5/scripts/add_electricity.py#L221-L259>`__. Demand can be chosen for weather years and socioeconomic `ssp` scenarios

* Network is built, cleaned and processed solely on `OpenStreetMap data <https://github.com/pypsa-meets-earth/pypsa-earth/blob/9acf89b8756bb60d61460c1dad54625f6a67ddd5/scripts/osm_pbf_power_data_extractor.py>`__

* Voronoi regions, where data is aggregated towards, can be replaced by administrative `GADM zones <https://github.com/pypsa-meets-earth/pypsa-earth/commit/4aa21a29b08c4794c5e15d4209389749775a5a52>`__

* `Augmented line expansion feature <https://github.com/pypsa-meets-earth/pypsa-earth/pull/175>`__ can make network meshed, connect isolated mini-grids to the main-grid.

* Community moved to `Discord <https://discord.gg/AnuJBk23FU>`__.

* Most meeting and agenda's are `open <https://github.com/pypsa-meets-earth/pypsa-earth#get-involved>`__.


Release Process
===============

* Checkout a new release branch ``git checkout -b release-v0.x.x``.

* Finalise release notes at ``doc/release_notes.rst``.

* Update ``envs/environment.fixed.yaml`` via
  ``conda env export -n pypsa-earth -f envs/environment.fixed.yaml --no-builds``
  from an up-to-date `pypsa-earth` environment.

* Update version number in ``doc/conf.py`` and ``*config.*.yaml``.

* Open, review and merge pull request for branch ``release-v0.x.x``.
  Make sure to close issues and PRs or the release milestone with it (e.g. closes #X).

* Tag a release on Github via ``git tag v0.x.x``, ``git push``, ``git push --tags``. Include release notes in the tag message.

* Upload code to `zenodo code repository <https://doi.org>`_ with `GPLv3 license <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.

* Create pre-built networks for ``config.default.yaml`` by running ``snakemake -j 1 extra_components_all_networks``.

* Upload pre-built networks to `zenodo data repository <https://doi.org/10.5281/zenodo.3601881>`_ with `CC BY 4.0 <https://creativecommons.org/licenses/by/4.0/>`_ license.

* Send announcement on the `PyPSA-Earth Discord channel <https://discord.gg/AnuJBk23FU>`_.
