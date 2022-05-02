..
  SPDX-FileCopyrightText: 2021 The PyPSA-Africa Authors

  SPDX-License-Identifier: CC-BY-4.0

##########################################
Release Notes
##########################################


Upcoming Release
================

**New Features and major Changes**

* Bug fixing (script retrieve_databundle) and rule run_test to ease testing `PR #322 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/322>`__

* Handling non-numerical entries in raw OSM data: `PR #287 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/287>`__

* General user experience improvements: `PR #326 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/326>`__

PyPSA-Africa 0.0.2 (6th April 2022)
=====================================

**New Features and major Changes**

* Plotting and summary features: `PR #211 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/211>`__ and `PR #214 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/214>`__

* Templates for issue, PR, feature request: `PR #216 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/216>`__

* Attach hydro enabled with all hydro types: `PR #232 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/232>`__

* Parallel download of osm data: `PR #232 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/232>`__

* Decoupling iso coding from geofabrik; rule download_osm_data extended to the world: `PR #236 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/236>`__

* Rule build_shape extended to the world: `PR #236 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/236>`__

* Validation of geofabrik links: `PR #249 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/249>`__

* Generalized version of Data retrieval with google and zenodo hosting platforms: `PR #242 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/242>`__ and `PR #260 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/260>`__

* Fix random state for kmean clustering, adopted from `PR 313 <https://github.com/PyPSA/pypsa-eur/pull/313>`__

* Implement area exclusions based on land type using the Copernicus Land Cover: `PR #272 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/272>`__.

* Flexible demand extraction for multiple years across the globe: `PR #275 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/275>`_

* Add CI caching and windows CI: `Commit CI windows <https://github.com/pypsa-meets-africa/pypsa-africa/commit/c98cb30e828cfda17692b8f5e1dd8e39d33766ad>`__,  `PR #277 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/277>`__.

* Change config to allow weather year extraction from snapshots as default: `PR #301 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/301>`__.

* Replace Restyler by .pre-commit `PR #307 https://github.com/pypsa-meets-africa/pypsa-africa/pull/307`__.

* Solved the issue of "overpassing nodes" and restyling osm_build_network: `PR #294 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/294>`__

* Revise deprecations in build_shape: `PR #315 <https://github.com/pypsa-meets-africa/pypsa-africa/pull/315>`__

PyPSA-Africa 0.0.1 (24th December 2021)
=====================================

This is the first release of PyPSA-Africa which heavily builds on `PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`__.

**New Features and major Changes**

* Include new data streams for Africa model

* Demand data implementation from `GEGIS <https://github.com/pypsa-meets-africa/pypsa-africa/blob/9acf89b8756bb60d61460c1dad54625f6a67ddd5/scripts/add_electricity.py#L221-L259>`__. Demand can be chosen for weather years and socioeconomic `ssp` scenarios

* Network is built, cleaned and processed solely on `OpenStreetMap data <https://github.com/pypsa-meets-africa/pypsa-africa/blob/9acf89b8756bb60d61460c1dad54625f6a67ddd5/scripts/osm_pbf_power_data_extractor.py>`__

* Voronoi regions, where data is aggregated towards, can be replaced by administrative `GADM zones <https://github.com/pypsa-meets-africa/pypsa-africa/commit/4aa21a29b08c4794c5e15d4209389749775a5a52>`__

* `Augmented line expansion feature <https://github.com/pypsa-meets-africa/pypsa-africa/pull/175>`__ can make network meshed, connect isolated mini-grids to the main-grid.

* Community moved to `Discord <https://discord.gg/AnuJBk23FU>`__.

* Most meeting and agenda's are `open <https://github.com/pypsa-meets-africa/pypsa-africa#get-involved>`__.


Release Process
===============

* Checkout a new release branch ``git checkout -b release-v0.x.x``.

* Finalise release notes at ``doc/release_notes.rst``.

* Update ``envs/environment.fixed.yaml`` via
  ``conda env export -n pypsa-eur -f envs/environment.fixed.yaml --no-builds``
  from an up-to-date `pypsa-eur` environment.

* Update version number in ``doc/conf.py`` and ``*config.*.yaml``.

* Open, review and merge pull request for branch ``release-v0.x.x``.
  Make sure to close issues and PRs or the release milestone with it (e.g. closes #X).

* Tag a release on Github via ``git tag v0.x.x``, ``git push``, ``git push --tags``. Include release notes in the tag message.

* Upload code to `zenodo code repository <https://doi.org>`_ with `GPLv3 license <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.

* Create pre-built networks for ``config.default.yaml`` by running ``snakemake -j 1 extra_components_all_networks``.

* Upload pre-built networks to `zenodo data repository <https://doi.org/10.5281/zenodo.3601881>`_ with `CC BY 4.0 <https://creativecommons.org/licenses/by/4.0/>`_ license.

* Send announcement on the `PyPSA-Africa Discord channel <https://discord.gg/AnuJBk23FU>`_.
