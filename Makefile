# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

.PHONY: test setup clean

test:
	set -e
	snakemake solve_all_networks -call --configfile config.tutorial.yaml # this runs the tutorial config
	snakemake solve_all_networks -call --configfile config.tutorial.yaml test/config.custom.yaml # add custom config to tutorial config
	snakemake solve_all_networks -call --configfile config.tutorial.yaml configs/scenarios/config.NG.yaml
	snakemake solve_all_networks_monte -call --configfile config.tutorial.yaml test/config.monte_carlo.yaml
	snakemake solve_all_networks -call --configfile config.tutorial.yaml test/config.landlock.yaml
	snakemake -c4 solve_sector_networks --configfile config.tutorial.yaml test/config.sector.yaml
	snakemake -c4 solve_sector_networks_myopic --configfile config.tutorial.yaml test/config.myopic.yaml
	echo "All tests completed successfully."

setup:
	# Add setup commands here
	echo "Setup complete."

clean:
	# Add clean-up commands here
	snakemake -j1 solve_all_networks --delete-all-output --configfile config.tutorial.yaml test/config.custom.yaml
	snakemake -j1 solve_all_networks --delete-all-output --configfile config.tutorial.yaml configs/scenarios/config.NG.yaml
	snakemake -j1 solve_all_networks_monte --delete-all-output --configfile test/config.monte_carlo.yaml
	snakemake -j1 run_all_scenarios --delete-all-output --configfile test/config.landlock.yaml
	snakemake -j1 solve_sector_networks --delete-all-output --configfile test/config.sector.yaml
	snakemake -j1 solve_sector_networks_myopic --delete-all-output --configfile config.tutorial.yaml test/config.myopic.yaml
	echo "Clean-up complete."
