<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->
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