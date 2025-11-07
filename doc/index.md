# Welcome to the PyPSA-Earth documentation!

[![GitHub release](https://img.shields.io/github/v/release/pypsa-meets-earth/pypsa-earth?include_prereleases)](https://github.com/pypsa-meets-earth/pypsa-earth/releases)
[![CI](https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/test.yml/badge.svg)](https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/test.yml)
[![Documentation Status](https://readthedocs.org/projects/pypsa-earth/badge/?version=latest)](https://pypsa-earth.readthedocs.io/en/latest/?badge=latest)
[![GitHub repo size](https://img.shields.io/github/repo-size/pypsa-meets-earth/pypsa-earth)](https://github.com/pypsa-meets-earth/pypsa-earth)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPLv3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![REUSE status](https://api.reuse.software/badge/github.com/pypsa/pypsa-eur)](https://api.reuse.software/info/github.com/pypsa/pypsa-eur)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/pypsa-meets-earth/pypsa-earth/main.svg)](https://results.pre-commit.ci/latest/github/pypsa-meets-earth/pypsa-earth/main)
[![Discord](https://img.shields.io/discord/911692131440148490?logo=discord)](https://discord.gg/AnuJBk23FU)
[![Google Drive](https://img.shields.io/badge/Google%20Drive-4285F4?style=flat&logo=googledrive&logoColor=white)](https://drive.google.com/drive/folders/1U7fgktbxlaGzWxT2C0-Xv-_ffWCxAKZz)

## Motivation

Closed-source models are the current standard for most policy and industry decisions. However, open models have proven to be competitive alternatives that promote science, robust technical analysis, collaboration and transparent policy decision making. Yet, two issues slow the adoption: open models are often designed with limited geographic scope, hindering synergies to collaborate, or are based on low spatially resolved data, limiting their utility.

PyPSA-Earth is the first open-source global cross-sectoral energy system model with high spatial and temporal resolution. The workflow provides capabilities for modelling the energy systems of any country in the world, enabling large-scale collaboration and transparent analysis for an inclusive and sustainable energy future. PyPSA-Earth is suitable for both operational studies and capacity expansion studies. Its sector-coupled modeling capabilities enable features for the detailed optimization of multi-energy systems, covering electricity, heating, transport, industry, hydrogen and more.

The **PyPSA meets Earth initiative** members are maintaining the *PyPSA-Earth* repository as well as many other tools. The [website](https://pypsa-meets-earth.github.io/) provides more context of the initiative and the associated projects.

<div class="grid-images" style="display: flex; gap: 3px; justify-content: center; align-items: stretch;">
  <img src="https://forum.openmod.org/uploads/db8804/original/1X/ddf041d1b98ca8f8c310f1c6393ec426ab5594cf.png" alt="Example 1" style="width: 32%; height: 300px;">
  <img src="https://forum.openmod.org/uploads/db8804/original/1X/940b2673cfc31c4a6f01b7908f546d39d67df27e.png" alt="Example 2" style="width: 32%; height: 300px;">
  <img src="https://forum.openmod.org/uploads/db8804/original/1X/6af089c376b19b72ad148e4e4326c162b94db68f.png" alt="Example 3" style="width: 32%; height: 300px;">
</div>

**Figure:** Example power systems built with PyPSA-Earth. See images of ~193 more countries at [https://zenodo.org/records/10080766](https://zenodo.org/records/10080766)

## Get Involved

There are multiple ways to get involved and learn more about our work:

1. **Join our forum** and communication platform on [PyPSA-meets-Earth Discord Server](https://discord.gg/AnuJBk23FU)

2. **Chat on Discord with us** in the following meetings:
    - General initiative meeting for project news and [high-level code updates](https://docs.google.com/document/d/1r6wm2RBe0DWFngmItpFfSFHA-CnUmVcVTkIKmthdW3g/edit?usp=sharing). Held every [fourth Thursday 16-17:00 (UK time)](https://drive.google.com/file/d/1naH4WwW9drkOkOJ3PLO4fyWdkZQi5-_w/view?usp=share_link) and is a perfect place to meet the community and get a high-level update on PyPSA ecosystem relevant for PyPSA-Earth developments.
    - Weekly developers meetings:
        - Eastern-Hemisphere friendly *Morning meeting* every [Thursday at 09:00 (UK time)](https://drive.google.com/file/d/1PDdmjsKhzyGRo0_YrP4wPQkn2XTNh6jA/view?usp=share_link).
        - Western-Hemisphere friendly *Evening meeting* every [Thursday 16:00 (UK time)](https://drive.google.com/file/d/1gaLmyV4qGPXsogkeRcAPWjC0ESebUxU-/view?usp=share_link). Every fourth Thursday is replaced by the General initiative meeting which has a more high-level perspective, but you can also join to discuss more particular questions.

3. **Look at public materials** at [Google Drive](https://drive.google.com/drive/folders/13Z8Y9zgsh5IZaDNkkRyo1wkoMgbdUxT5?usp=sharing) to share minutes, presentations, lists and documents. Feel free to get a look!

4. **Notify your interest** to on-demand meetings:
    - Demand creation and prediction meeting
    - AI asset detection meeting
    - Outreach meeting for planning, discussing events, workshops, communication, community activities

5. Join us and **propose your stream**.

## Documentation

### Getting Started

* [Introduction](home/introduction.md) - What is PyPSA-Earth and why use it
* [Installation](home/installation.md) - Set up your environment
* [Quick Start](home/quick-start.md) - Run your first model

### Tutorials

* [Electricity Model](tutorials/electricity-model.md) - Build an electricity-only model
* [Sector-Coupled Model](tutorials/sector-coupled-model.md) - Create a multi-sector model
* [Examples](tutorials/examples.md) - Jupyter notebooks and use cases

### User Guide

* [Configuration](user-guide/configuration.md) - Configure your model settings
* [Wildcards](user-guide/wildcards.md) - Understand wildcard patterns
* [Costs](user-guide/costs.md) - Technology cost assumptions
* [Structure](user-guide/structure.md) - Project structure and workflow
* [Rules Overview](user-guide/rules-overview.md) - Snakemake rules explained
* [Optimization](user-guide/optimization.md) - Optimization theory and methods
* [Plotting](user-guide/plotting.md) - Visualization and results
* [Data Workflow](user-guide/data-workflow.md) - Data processing pipeline
* [Model Customization](user-guide/model-customization.md) - Customize your model

### Advanced

* [Monte Carlo](advanced/monte-carlo.md) - Uncertainty quantification and sensitivity analysis

### API Reference

* [API Documentation](api-reference/index.md) - Complete API reference

### Community & Resources

* [Contributing](community/contributing.md) - How to contribute to the project
* [Project Structure](community/project-structure.md) - Credits and architecture
* [Release Notes](community/release-notes.md) - Version history and changes
* [Users List](community/users-list.md) - Who's using PyPSA-Earth
* [Talks & Papers](community/talks-and-papers.md) - Publications and presentations
* [Learning Materials](community/learning-materials.md) - Additional resources
* [Software Hints](community/software-hints.md) - Tips and troubleshooting
