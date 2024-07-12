# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

.PHONY: test setup clean

test:
	set -e
	snakemake solve_all_networks -call # this runs the tutorial config
	snakemake solve_all_networks -call --configfile test/config.custom.yaml # add custom config to tutorial config
	snakemake solve_all_networks -call --configfile configs/scenarios/config.NG.yaml
	snakemake solve_all_networks_monte -call --configfile test/config.monte_carlo.yaml
	snakemake solve_all_networks -call --configfile test/config.landlock.yaml
	echo "All tests completed successfully."

setup:
	# Add setup commands here
	echo "Setup complete."

clean:
	# Add clean-up commands here
	snakemake -j1 solve_all_networks --delete-all-output --configfile test/config.custom.yaml
	snakemake -j1 solve_all_networks --delete-all-output --configfile test/config.tutorial_noprogress.yaml
	snakemake -j1 solve_all_networks_monte --delete-all-output --configfile test/config.monte_carlo.yaml
	snakemake -j1 run_all_scenarios --delete-all-output --configfile test/config.landlock.yaml
	echo "Clean-up complete."
