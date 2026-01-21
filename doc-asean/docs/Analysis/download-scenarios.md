# Download Scenarios

The PyPSA scenario results and sensitivity models can be downloaded here:

- [results-asean-paper-v1.zip](https://drive.google.com/file/d/194my_d3eotQ1GvsGGN5KOZGWOHEivoQy/view?usp=drive_link)

There are multiple options for exploring and analyzing the model results.

## Jupter Notebooks

After downloading the data, ensure that the files are placed in the following directory structure:

```
results/
├── baseline-aims-3H/
│   └── postnetworks/
│       ├── elec_s_100_ec_lv2.0__3h_2025_0.071_DEC_0export.nc
│       ├── ...
│       └── elec_s_100_ec_lv2.0__3h_2050_0.071_DEC_0export.nc
├── ...
└── sensitivity-roof200-3H/

```

Descriptions and usage instructions for each notebook are provided [here](jupyter-notebooks.md)

## PyPSA App

[PyPSA App](https://github.com/PyPSA/pypsa-app) is a self-hosted web application for analyzing and visualizing PyPSA networks, with modular architecture designed to extend into workflow execution, network editing, optimization, and custom integrations.

> **Warning:** This app is in early development.

