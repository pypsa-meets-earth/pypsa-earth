<!--
SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
SPDX-License-Identifier: CC-BY-4.0
-->

# Alternative Clustering

Alternative clustering is not a standalone Snakemake rule but a configuration option within the `cluster_options` section that changes the behavior of several rules, including `cluster_network`, `build_renewable_profiles`, `build_powerplants`, and `add_electricity`.

## What does it do?

When `alternative_clustering` is enabled, the model uses administrative boundaries (GADM shapes) instead of Voronoi regions for spatial clustering. This can be useful when:

- The default Voronoi-based clustering does not align well with administrative or political boundaries
- You want to model specific regions or provinces as individual network nodes
- Regional policy analysis requires alignment with administrative divisions

## Configuration

```yaml
cluster_options:
  alternative_clustering: false  # Set to true to enable GADM-based clustering
```

## Affected Rules

The `alternative_clustering` flag is passed as a parameter to several rules:

| Rule | Effect |
|---|---|
| `build_renewable_profiles` | Changes how renewable profiles are assigned to regions |
| `build_powerplants` | Changes how power plants are assigned to buses |
| `add_electricity` | Adjusts how generators are attached to the network |
| `cluster_network` | Uses GADM shapes for clustering instead of Voronoi regions |
| `prepare_gas_network` | Adjusts gas network clustering |
| `prepare_sector_network` | Adjusts sector network clustering |
