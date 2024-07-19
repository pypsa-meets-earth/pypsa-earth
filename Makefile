# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

.PHONY: test setup clean

test:
	set -e
	snakemake -call solve_sector_networks --configfile pypsa-earth/config.tutorial.yaml test/config.test1.yaml
	echo "All tests completed successfully."

setup:
	# Add setup commands here
	echo "Setup complete."

clean:
	# Add clean-up commands here
	snakemake -j1 solve_sector_networks --delete-all-output --configfile test/config.test1.yaml
	echo "Clean-up complete."
