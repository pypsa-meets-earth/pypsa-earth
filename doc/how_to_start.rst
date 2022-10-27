..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _how_to_start:


##########################################
How to run your model
##########################################

The installation procedure installs PyPSA-Earth model with all the dependencies needed to build it. The next step is to introduce the data needed to model a region of your interest. Flexible data management is the core feature of PyPSA-Earth model which relies on customizable data extraction and preparation scripts with global coverage.

Tutorial run
------------------------------

The entire distribution, including the data for most parts on Earth, is very heavy (e.g. for Africa it is >40Gb) and it involves a large number of files. To facilitate model testing, a lightweight data starter kit was developed. You can use it by setting `tutorial: true` `retrieve_databundle: true` in the configuration file `config.yaml` placed in the project folder `pypsa-earth`. 

After doing so, you just need run the model with the following command:

.. code:: bash

    snakemake -j 1 solve_all_networks

.. TODO Explain settings of the tutorial case

Working run
------------------------------

After playing with the tutorial model, it's important to clean-up data in your model folder to avoid data conflicts. You may use the `clean` rule for making so:

.. code:: bash

    snakemake -j 1 clean

Generally, it's a good idea to repeat the cleaning procedure every time when the underlying data are changed.

It's recommended to set `retrieve_databundle: true` when building the model first time to download the common data files needed. The load will start automatically when running the model with:

.. code:: bash
    snakemake -j 1 solve_all_networks

Please, note that around **20 Gb zipped files will be downloaded**. It's worth to make sure you have a stable connection, time and around 50 Gb available in your system. If no errors show up, then you can proceed. It's advicable to set `retrieve_databundle: false` after the first model run when all the needed data will be successfully extracted to avoid data loss.

Snakemake
===========================

Snakemake is a workflow management tool inherited by PyPSA-Earth from PyPSA-Eur. Snakemake decomposes a large software process into a set of subtasks, or ’rules’, that are automatically chained to obtain the desired output.

.. note::
  ``Snakemake``, which is one of the major dependencies, will be automatically installed in the environment pypsa-earth, thereby there is no need to install it manually.

The snakemake included in the conda environment pypsa-earth installed with the above-mentioned procedure can be executed with the following procedure:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake < any command here >  

Starting with essential usability features, the implemented Snakemake procedure allows to flexibly execute the entire workflow with various options without writing a single line of python code. For instance, you can model the world energy system or any subset of countries only using the required data. Wildcards, which are special generic keys that can assume multiple values depending on the configuration options, help to execute large workflows with parameter sweeps and various options.


.. TODO Add Snakemake tutorial links    
