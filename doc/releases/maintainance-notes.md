<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->
# Release Process

* Make sure to merge all PRs intended to be included into a new PR and communicate to the authors of the linked issues that the respective PR have been completed

* Checkout the `main` branch locally and pull the recent changes locally from `PyPSA-Earth` upstream.

* Checkout a new release branch `git checkout -b release-v0.x.x`.

* Review and finalise release notes at `doc/release_notes.rst`.

* Create a draft release using GitHub interface with `Draft release` under https://github.com/pypsa-meets-earth/pypsa-earth/releases. Review and fix the release notes.

* Make sure that locked versions of the environments `*-plock.yaml` in `envs` folder, as well as `pixi.lock` are up-to-date.

* Update version number in `default.config.yaml`, `tutorial.config.yaml` and `test/config.*.yaml`.

* Open, review and merge pull request for branch `release-v0.x.x`. Run `pre-commit run --all`` locally and fix any issues.

* Update and checkout your local `main` and tag a release with `git tag v0.x.x`, `git push`, `git push --tags`. Include release notes in the tag message using Github UI.

* Send announcement on the [PyPSA-Earth Discord channel](https://discord.gg/AnuJBk23FU) and LinkedIn.
