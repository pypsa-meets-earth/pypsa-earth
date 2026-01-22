# Jupyter Notebooks

This document briefly describes the Jupyter notebooks in the `notebooks/` folder and gives quick usage notes so you can reproduce the figures and tables used in the PyPSA-ASEAN analysis.

To use this notebook, you need a PyPSA version of at least 1.0. To create a separate conda environment from an existing one:

```bash
$ conda create --name pypsa-earth-lab --clone pypsa-earth
$ conda activate pypsa-earth-lab
$ conda update -c conda-forge pypsa>=1.0
```

Then, activate jupyter lab:
```bash
$ jupyter lab
```

The `analysis_helper.py` contains dictionary and functions that are used across all notebooks mentioned below. To save the figure results, set `savefig = True` in the respective notebooks.

## Preparation Analysis

The `preparation_analysis.ipynb` generates plots of networks prior to the solving process. This includes visualizing renewable energy potentials and displaying existing energy network in the region.

## Validation Analysis

The `validation_analysis.ipynb` compares the baseline scenario for the year 2025 with historical results based on EMBER. Before using this script, you must have at least one network that has been solved beforehand.

## Postnetwork Analysis

The `postnetwork_analysis.ipynb` generates an overview of the network after optimization. These scripts use features from PyPSA v1.0 to create KPI tables, energy network maps, and figures showing the total energy mix and system cost.

## Regional Analysis

The `regional_analysis.ipynb` focuses on the regional dimension of the scenario results by disaggregating the energy mix by country and subregion. Additionally, there is an unused script for regional marginal prices.

## Sensitivity Analysis

The `sensitivity_analysis.ipynb` compares the results of sensitivity scenarios with those of a selected benchmark scenario in terms of total energy mix and system cost.

