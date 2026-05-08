<!--
SPDX-FileCopyrightText:  PyPSA-Earth, PyPSA-Zambia and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

# Welcome to the PyPSA-Zambia documentation

[![GitHub release](https://img.shields.io/github/v/release/open-energy-transition/pypsa-zambia?include_prereleases)](https://github.com/open-energy-transition/pypsa-zambia/releases)
[![CI](https://github.com/open-energy-transition/pypsa-zambia/actions/workflows/test.yml/badge.svg)](https://github.com/open-energy-transition/pypsa-zambia/actions/workflows/test.yml)
[![Documentation Status](https://readthedocs.org/projects/pypsa-zambia/badge/?version=latest)](https://pypsa-zambia.readthedocs.io/en/latest/?badge=latest)
[![GitHub repo size](https://img.shields.io/github/repo-size/open-energy-transition/pypsa-zambia)](https://github.com/open-energy-transition/pypsa-zambia)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPLv3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![REUSE status](https://api.reuse.software/badge/github.com/open-energy-transition/pypsa-zambia)](https://api.reuse.software/info/github.com/open-energy-transition/pypsa-zambia)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/open-energy-transition/pypsa-zambia/main.svg)](https://results.pre-commit.ci/latest/github/open-energy-transition/pypsa-zambia/main)
[![Discord](https://img.shields.io/discord/911692131440148490?logo=discord)](https://discord.gg/AnuJBk23FU)
[![Google Drive](https://img.shields.io/badge/Google%20Drive-4285F4?style=flat&logo=googledrive&logoColor=white)](https://drive.google.com/drive/folders/1U7fgktbxlaGzWxT2C0-Xv-_ffWCxAKZz)

## Open Tools to Improve Grid Planning and Operation in Africa

Zambia’s national electricity utility ZESCO and the non-profit Open Energy Transition are developing an open-source energy modelling tool
aimed at improving data-driven power system planning in Zambia. The project will support the national grid operators and utilities in adopting transparent, climate-resilient planning approaches while strengthening local technical capacity.

Open Energy Transition are developing this customized energy system model based on PyPSA (Python for Power System Analysis), tailored to the Zambian energy system and designed with scaling potential for other African utilities. Tool development will be guided by local insights to ensure the resulting solutions are practical, accessible, and trusted by those who rely on them. Other partners include Association of Power Utilities in Africa (APUA), the Zambian Ministry of Energy, and academic partners.

The collaboration will run for 36 months. Initial concept development began in Q4 2025 and is now in implementation.
The project is funded by the Quadrature Climate Foundation.

As Africa’s electricity systems expand rapidly, planning decisions made today will shape reliability, affordability, and emissions for decades. Zambia is an ideal pilot country as energy shocks and multi-day blackouts have highlighted the urgent need for resilient energy planning. Meanwhile, the African Continental Master Plan depends on local capacity to develop national plans.

The open-source modelling framework PyPSA is particularly well suited to this context. It is cost-efficient, transparent, and reproducible – ideal to tackle the complexity of modern energy systems and the move towards a just and sustainable energy transition.


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

- [Introduction](home/introduction.md) - What is PyPSA-Earth and why use it
- [Installation](home/installation.md) - Set up your environment
- [Quick Start](home/quick-start.md) - Run your first model

### Tutorials

- [Electricity Model](tutorials/electricity-model.md) - Build an electricity-only model
- [Sector-Coupled Model](tutorials/sector-coupled-model.md) - Create a multi-sector model
- [Model Customization](user-guide/model-customization.md) - General Modeling Guidelines
- [Examples](tutorials/examples.md) - Jupyter notebooks and use cases

### User Guide

- [Structure](user-guide/structure.md) - Project structure and workflow
- [Data Workflow](user-guide/data-workflow.md) - Data processing pipeline
- [Configuration](user-guide/configuration.md) - Configure your model settings
- [Wildcards](user-guide/wildcards.md) - Understand wildcard patterns
- [Costs](user-guide/costs.md) - Technology cost assumptions
- [Rules Overview](user-guide/rules-overview.md) - Snakemake rules explained
- [Rules Reference](user-guide/rules-reference/download-and-filter/index.md) - Detailed rule descriptions

### Educational Materials

- [Optimization](user-guide/optimization.md) - Optimization theory and methods

### Utilities

- [Monte Carlo](utilities/monte-carlo.md) - Uncertainty quantification and sensitivity analysis

### Community & Resources

- [Contributing](community/contributing.md) - How to contribute to the project
- [Project Structure](community/project-structure.md) - Credits and architecture
- [Users List](community/users-list.md) - Who's using PyPSA-Earth
- [Talks & Papers](community/talks-and-papers.md) - Publications and presentations
- [Learning Materials](community/learning-materials.md) - Additional resources
- [Software Hints](community/software-hints.md) - Tips and troubleshooting

### API Reference

- [API Documentation](api-reference/index.md) - Complete API reference
