<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->
# Release Process

* Checkout the `main` branch locally and pull the recent changes locally from `PyPSA-Earth` upstream.

* Checkout a new release branch `git checkout -b release-v0.x.x`.

* Review and finalise release notes at `doc/release_notes.rst`.

* Make sure that locked versions of the environments `*-plock.yaml` in `envs` folder, as well as `pixi.lock` are up-to-date.

* Update version number in `default.config.yaml`, `tutorial.config.yaml` and `test/config.*.yaml`.

* Open, review and merge pull request for branch `release-v0.x.x`.
  Make sure to close issues and PRs or the release milestone with it (e.g. closes #X).
  Run `pre-commit run --all`` locally and fix any issues.

* Update and checkout your local `main` and tag a release with `git tag v0.x.x`, `git push`, `git push --tags`. Include release notes in the tag message using Github UI.

* Upload code to `zenodo code repository](<https://doi.org>) with [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html).

* Send announcement on the [PyPSA-Earth Discord channel](https://discord.gg/AnuJBk23FU).
