# Dual Pull Request

This section details how PyPSA-ASEAN employs a "dual pull request" merge strategy to efficiently integrate globally relevant features. This approach facilitates rapid implementation within PyPSA-ASEAN while ensuring that features are thoroughly reviewed and eventually merged into the upstream PyPSA-Earth project.

## How it Works

The strategy involves the following steps:

1. A gap in required features for PyPSA-ASEAN is identified.
2. The identified feature is determined to have broader applicability beyond the ASEAN region.
3. A pull request is simultaneously made from a fork of PyPSA-Earth to both the PyPSA-ASEAN repository and the upstream PyPSA-Earth repository.
4. If the feature functions as intended, it is merged into PyPSA-ASEAN for immediate use.
5. Further improvements and adjustments are made based on suggestions and reviews from the PyPSA-Earth community.
6. Once refined and approved, the pull request is merged into PyPSA-Earth.
7. PyPSA-ASEAN, following a soft-fork strategy, continuously updates and benefits from the latest features integrated into PyPSA-Earth.

## Implemented Features

The following features have been implemented using this dual-pull-request method:

*   Add psuedo `branch()` to streamline snakemake workflow ([PR #10](https://github.com/pypsa-meets-earth/pypsa-asean/pull/10)).
*   Add subregionalization in `cluster_networks` ([PR #11](https://github.com/pypsa-meets-earth/pypsa-asean/pull/11)).
*   Add `co2_budget` for emission targets in multiple planning horizon years ([PR #13](https://github.com/pypsa-meets-earth/pypsa-asean/pull/13)).
*   Reintroduce `sanitize_carriers` and `sanitize_location` ([PR #14](https://github.com/pypsa-meets-earth/pypsa-asean/pull/14)).
*   Use Global Buildings datasets to estimate solar-rooftop potentials ([PR #20](https://github.com/pypsa-meets-earth/pypsa-asean/pull/20)).
*   Bugfix: Avoid creating duplicate conventional generators by setting k`eep_existing_capacities` to `false` ([PR #21](https://github.com/pypsa-meets-earth/pypsa-asean/pull/21)).
*   Bugfix: fix missing baseyear capacity, set biomass as default for bioenergy ([PR #22](https://github.com/pypsa-meets-earth/pypsa-asean/pull/22)).
*   Add IRENA statistics for hydro generation normalization ([PR #23](https://github.com/pypsa-meets-earth/pypsa-asean/pull/23)).
*   Add heuristics to determine missing hydro technologies ([PR #24](https://github.com/pypsa-meets-earth/pypsa-asean/pull/24)).