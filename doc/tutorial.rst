..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Africa authors

  SPDX-License-Identifier: CC-BY-4.0

.. _tutorial:

##########################################
Tutorial
##########################################

.. _prerequisites_learning_material:

Prerequisites and learning material
===================================

PyPSA meets Africa builds on top of several open-source packages, which are here recalled together with recommended sources to learn them from scratch.

.. _data_science_basics:

Data science basics (essential)
--------------------------------


- Refresh your Python knowledge by watching `CSDojo's playlist <https://www.youtube.com/c/CSDojo/playlists>`_. His content is excellent as introduction. You will learn in effective short videos the python basics such as variables If/else statements, functions, lists, for loops, while loops, dictionaries, classes and objects, boolean, list comprehensions, sets - put your hands on and write some test scripts as the video suggests. (~3h)
- Familiarize yourself with numpy and panda dataframes.  In the Python-based PyPSA tool, we do not work with Excel. Powerful panda dataframes are our friends. `Here <https://www.coursera.org/learn/python-data-analysis>`__ is an extensive 30h course that provides a great introduction if this is completely unfamiliar to you.
- `Introduction to Unix-shell <https://swcarpentry.github.io/shell-novice/>`_ - "Use of the shell is fundamental to a wide range of advanced computing tasks, including high-performance computing and automated workflow. These lessons will introduce you to this powerful tool." (optional 4h, to become a pro)


PyPSA Introduction (essential)
-------------------------------

- Watch how PyPSA-Eur is designed https://www.youtube.com/watch?v=ty47YU1_eeQ (1h)
- Watch and put your hands on to make PyPSA-Eur work on your computer https://www.youtube.com/watch?v=mAwhQnNRIvs (1-3h)
- While watching these PyPSA videos always have a look into the excellent `PyPSA-Eur documentation <https://pypsa-eur.readthedocs.io/en/latest/index.html>`_
- To see what data we can extract we work usually closely with the `basic PyPSA documentation <https://pypsa.readthedocs.io/en/latest/components.html>`_ 


Git and GitHub (essential)
---------------------------

For code collaboration we use GitHub. Which is a common source control tool that is a very popular collaborative code development tool. Here some notes if you are not already familar with it:

- Git and GitHub is not the same. Usually, you work with git on your computer (offline) to push changes to GitHub (online).
- `Here <https://www.youtube.com/watch?v=8JJ101D3knE>`__ a great intro which we recommend
- Learning by doing. Maybe one of the best ways to learn is to puts your hands on open a GitHub repository and upload/change/reverse files from your local computer on some dummy scripts.
- This `cheatsheet <https://www.atlassian.com/git/tutorials/atlassian-git-cheatsheet>`_ might help using the Git commands


Snakemake and advanced changes (essential)
-------------------------------------------

Snakemake is our brain in PyPSA. 
It automates many tasks & keeps the code structure clean. 
Therefore, it is quite useful to learn if your task is to integrate features into PyPSA.
We can recommend: 

- `snakemake basic and advanced tutorial here <https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html>`__ (takes max 3-5h and makes a lot of fun).
- Explore how PyPSA uses snakemake in the Snakefile and scripts - the GitHub search function is your best friend to find quickly what interests you.


Code environment (optional)
-----------------------------

We can recommend setting up VScode from Microsoft. Add some extension if you like as described in `this video <https://www.youtube.com/watch?v=0fROnrISdZU>`_. For instance GitHub, Gitlense, and maybe some others.

*Note*: if you decide to use Visual Studio Code, check out the tutorial about how to use `Git <https://code.visualstudio.com/docs/editor/versioncontrol#_git-support>`_ and `Github <https://code.visualstudio.com/docs/editor/github>`_  in Visual Studio Code


Notebooks
===========

In order to familiarize with the code or investigate the input and outputs of
the rules, in the ``notebooks`` folder, the following notebooks are available:

- ``network_comparison``: compares the network models developed along the data workflow; useful and interactive plots are generated
- ``osm_build_network_plot``: provides specific plots and outputs for the ``download_osm_data`` rule
- ``osm_data_access``: explains how OSM data are being loaded by using ``download_osm_data``
- ``osm_powermap``: contains nice plots and description of the output of the data downloaded and cleaned by using ``download_osm_data`` and ``clean_osm_data``
- ``solve_network_results``: provides useful plots and textual outputs to investigate the results of the last optimization performed using solve_network
- ``build_bus_regions``: it explores the inputs and outputs of the ``build_bus_regions`` rule,
  namely the bus regions shapes, and the elements of the network (lines, substations, etc.)
- ``build_shapes``: describes the shapes created using the rule ``build_shapes`` for the on-shore and off-shore areas
- ``demand_gegis``: it enables exploring th GeGIS dataset used to perform the analysis.
  These data are obtained using the `GlobalEnergyGIS <https://github.com/niclasmattsson/GlobalEnergyGIS>`_ package for Africa.
- ``shape_comparison``: this notebook enables comparing the shapes used along the data workflow
- ``add_electricity``: it analyzes the outputs of the ``add_electricity`` rule, including the PyPSA model and the RES/demand inputs
- ``base_network``: it eases the visualization and analysis of the output PyPSA network model that the rule ``base_network`` builds
- ``build_cutout``: the notebook analyzes the outputs of the rule ``build_cutout`` rule, which are the solar, wind and hydro time series
  generated with `Atlite <https://github.com/PyPSA/atlite/>`_
- ``build_renewable_profiles``: it enables investigating the specific time series generated by the rule ``build_renewable_profiles``;
  in particular, it shows the potential of selected resources (e.g. solar) and the corresponding time series of renewable energy production
  available for selected buses
- ``landuse-availability``: this notebook aims at showing how ``Atlite`` accounts for land constraints in the analysis


Examples
========

Solve the optimal power flow
-----------------------------------

The following snakemake routine enables executing the optimal power flow problem
for the current configuration using a 6-bus equivalent of the region.

.. code:: bash

    .../pypsa-africa % snakemake -j 1 results/networks/elec_s_6_ec_l.nc

Change the country for the analysis
-----------------------------------

In order to run the code for a set of countries different than the default ones,
the option ``countries`` in the configuration yaml files shall be modified.
To do so, follow the following procedure:

1. Make a copy of the ``config.default.yaml`` file and rename it as ``config.yaml``
2. In ``config.yaml`` modify the option ``countries = ["AA", ..., "ZZ"]`` with the list
   of countries that you desire; 2-digit country codes are requested or region names.
   
   For example, to investigate Nigeria, the following specification shall be applied in
   the configuration file.

   .. code:: bash

      countries = ["NG"]

   The code also supports pre-set group of countries, such as africa. For example,
   the African region can be simulated using:

   .. code:: bash

      countries = ["africa"]
3. Then, the software is ready to be used on the selected countries

Manual test of specific scripts
-------------------------------

The scripts in the ``scripts`` folder are build so that they can be easily run and tested
even without the snakemake procedure. Therefore, to test the specific functionality of
