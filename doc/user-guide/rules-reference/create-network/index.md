<!--
SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
SPDX-License-Identifier: CC-BY-4.0
-->

# Create Network

The simplification `snakemake` rules prepare **approximations** of the full model, for which it is computationally viable to co-optimize generation, storage and transmission capacities.

## Rules

- **[base_network](base-network.md)** - Builds and stores the base network with all buses, HVAC lines and HVDC links.

- **[add_electricity](add-electricity.md)** - Adds the generators and demand to the network model.

- **[simplify_network](simplify-network.md)** - Transforms the transmission grid to a 380 kV only equivalent network.

- **[cluster_network](cluster-network.md)** - Uses a clustering approach (e.g. [k-means](https://en.wikipedia.org/wiki/K-means_clustering)) to partition the network into a given number of zones and then reduce the network to a representation with one bus per zone.

- **[alternative_clustering](alternative-clustering.md)** - Provides alternative algorithms and methods for network clustering when the default algorithms are not optimal for the studied region.

- **[augmented_line_connections](augmented-line-connections.md)** - Reinforces or creates user-defined missing line connections to ensure the grid is fully connected and avoids isolated grid components.

- **[add_extra_components](add-extra-components.md)** - Adds extra components to the model, such as storage.

- **[prepare_network](prepare-network.md)** - Introduces optional constraints and requirements in the modelling, such as CO2 emissions, security margins, etc.

## Reference

The simplification and clustering steps are described in detail in the paper:

Jonas Hörsch and Tom Brown. [The role of spatial scale in joint optimisations of generation and transmission for European highly renewable scenarios](https://arxiv.org/abs/1705.07617), *14th International Conference on the European Energy Market*, 2017.
