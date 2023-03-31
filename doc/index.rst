.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. PyPSA meets Earth documentation master file, created by
   sphinx-quickstart on Sat May 15 22:52:54 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the PyPSA-Earth documentation!
================================================


.. image:: https://img.shields.io/github/v/release/pypsa-meets-earth/pypsa-earth?include_prereleases
    :alt: GitHub release (latest by date including pre-releases)

.. image:: https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/ci-linux.yaml/badge.svg
    :target: https://github.com/pypsa-meets-earth/pypsa-earth/actions
    :alt: CI Linux

.. image:: https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/ci-mac.yaml/badge.svg
    :target: https://github.com/pypsa-meets-earth/pypsa-earth/actions
    :alt: CI Mac

.. image:: https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/ci-windows.yaml/badge.svg
    :target: https://github.com/pypsa-meets-earth/pypsa-earth/actions
    :alt: CI Windows

.. image:: https://readthedocs.org/projects/pypsa-earth/badge/?version=latest
    :target: https://pypsa-earth.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://img.shields.io/github/repo-size/pypsa-meets-earth/pypsa-earth
    :alt: GitHub repo size

.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
    :target: https://www.gnu.org/licenses/gpl-3.0
    :alt: License

.. image:: https://api.reuse.software/badge/github.com/pypsa/pypsa-eur
    :target: https://api.reuse.software/info/github.com/pypsa/pypsa-eur
    :alt: REUSE

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: Code style Black

.. image:: https://results.pre-commit.ci/badge/github/pypsa-meets-earth/pypsa-earth/main.svg
    :target: https://results.pre-commit.ci/latest/github/pypsa-meets-earth/pypsa-earth/main
    :alt: Pre-commit CI-status

.. image:: https://img.shields.io/discord/911692131440148490?logo=discord
    :target: https://discord.gg/AnuJBk23FU
    :alt: Discord

.. image:: https://img.shields.io/badge/Google%20Drive-4285F4?style=flat&logo=googledrive&logoColor=white
    :target: https://drive.google.com/drive/folders/1U7fgktbxlaGzWxT2C0-Xv-_ffWCxAKZz
    :alt: Google Drive

*Motivation*. Closed-source models are the current standard for most policy and industry decisions. However, open models have proven to be
competitive alternatives that promote science, robust technical analysis, collaboration and transparent policy decision making.
Yet, two issues slow the adoption: open models are often designed with limited geographic scope, hindering synergies to collaborate,
or are based on low spatially resolved data, limiting their utility.

*PyPSA-Earth* is the first open-source global energy system model with data in high spatial and temporal resolution. It enables
large-scale collaboration by providing a tool that can model the world energy system or any subset of it. This work is derived
from the European PyPSA-Eur model using new data and functions. It is suitable for operational as well as combined generation,
storage and transmission expansion studies. We work hard to extend the PyPSA-Earth model by end of this year to include sector-coupling,
myopic and perfect pathway expansion capabilities.

*PyPSA meets Earth initiative* members are maintaining the *PyPSA-Earth* repository as well as many other tools.
The `website <https://pypsa-meets-earth.github.io/>`_ provides more context of the initiative and the associated projects. 

.. image:: img/africa_osm_map.png
    :width: 60%
    :align: center

==============
Get Involved
==============

The PyPSA meets Earth team is currently running four types of meetings:


- `**Discord NEW! (Open)** <https://discord.gg/AnuJBk23FU>`_

  - Chat with the community, team up on features, exchange with developers, code in voice channels

- **General code meeting (Open)**

  - every forth Thursday each month 16-17:00 (UK time) `download .ics file <https://drive.google.com/file/d/1naH4WwW9drkOkOJ3PLO4fyWdkZQi5-_w/view?usp=share_link>`_
  - updates on overall project and code blocks
  - meeting hosted on `Discord <https://discord.gg/AnuJBk23FU>`_; join us, we are waiting for you!
  - `open agenda <https://docs.google.com/document/d/1r6wm2RBe0DWFngmItpFfSFHA-CnUmVcVTkIKmthdW3g/edit?usp=sharing>`_. See what we will discuss. Invited members have edit rights.

- **Specific code meeting (by invitation)**

  - meeting hosted on Discord
  - join updates, demos, Q&A's, discussions and the coordination of each work package
  
    1. Demand creation and prediction meeting, on demand basis
    2. AI asset detection meeting, on demand basis
    3. Sector coupling meeting, every Thursday 09:00 (UK time), `download .ics file <https://drive.google.com/file/d/1PDdmjsKhzyGRo0_YrP4wPQkn2XTNh6jA/view?usp=share_link>`__
    4. PyPSA-Earth meeting, every Thursday 16:00 (UK time), `download .ics file <https://drive.google.com/file/d/1gaLmyV4qGPXsogkeRcAPWjC0ESebUxU-/view?usp=share_link>`__

- **Outreach meeting (by invitation)**

  - every second week
  - planning, discussing events, workshops, communication, community activities
  
- **Buddy talk (Open)**

  - book a 30min meeting with Max to discuss anything you like
  - booking link: `calendly.com/pypsa-meets-earth <https://calendly.com/pypsa-meets-earth/pypsa-meets-earth-exchange-30min>`_ (developed by @mnm-matin)
calendly.com/pypsa-meets-earth

=============
Documentation
=============

**Getting Started**

* :doc:`introduction`
* :doc:`installation`
* :doc:`short_tutorial`
* :doc:`tutorial`
* :doc:`data_workflow`
* :doc:`notebooks`

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Getting Started

   introduction
   installation  
   short_tutorial
   tutorial
   data_workflow    
   notebooks

**Configuration**

* :doc:`wildcards`
* :doc:`configuration`
* :doc:`costs`

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Configuration

   wildcards
   configuration
   costs

**Work flow and API**

* :doc:`structure`
* :doc:`rules_overview`
* :doc:`api_reference`

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Work flow and API

   structure
   rules_overview
   api_reference

**Help and References**

* :doc:`release_notes`
* :doc:`how_to_contribute`
* :doc:`software_hints`
* :doc:`learning_materials`
* :doc:`project_structure_and_credits`
* :doc:`talks_and_papers`

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Project Info

   release_notes
   how_to_contribute
   software_hints
   learning_materials
   project_structure_and_credits
   talks_and_papers
   
