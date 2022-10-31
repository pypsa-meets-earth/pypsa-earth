..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _quick_start:


##########################################
How to start with modeling
##########################################

The installation procedure installs PyPSA-Earth model with all the software dependencies needed to build and run it. To properly model any region of the Earth, PyPSA-Earth needs to download and fetch different datasets. This section explains how to perform this data management.

Build the tutorial model
------------------------------

The whole global data kit is very heavy (e.g. for Africa it is >40Gb) and involves a large number of files. To facilitate model testing, a lightweight tutorial kit was developed. You can use it by using the tutorial configuration file `config.tutorial.yaml` (placed in the project folder `pypsa-earth`). To do that, you may want to do a reserve copy of your current configuration file and then overwrite it by a tutorial configuration:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % cp config.tutorial.yaml config.yaml

cp config.tutorial.yaml config.yaml


In the configuration file `config.yaml` there is a flag `retrieve_databundle` which triggers data loading and a `tutorial` flag which determines that the loaded data belong to the light tutorial kit:

.. code:: yaml

    tutorial: true
    ...
    retrieve_databundle: true

How to run the model?
------------------------------

After configuration set-up, you just need run the modeling workflow with the following command:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 solve_all_networks

.. TODO Explain settings of the tutorial case

This command will trigger loading of the whole dataset needed to build the model for a tutorial case that is simulation of power systems in Nigeria and Benin. Note please that data loading and model building will take a while (about 20..50 minutes).

How to analyse the solved networks?
------------------------------

After playing with the tutorial model, it's important to clean-up data in your model folder to avoid data conflicts. You may use the `clean` rule for making so:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 clean

Generally, it's a good idea to repeat the cleaning procedure every time when the underlying data are changed.

It's recommended to set `retrieve_databundle: true` when building the model first time to download the common data files needed. The load will start automatically when running the model with:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 solve_all_networks

Please, note that around **20 Gb zipped files will be downloaded**. It's worth to make sure you have a stable connection, time and around 50 Gb available in your system. If no errors show up, then you can proceed. It's advisable to set `retrieve_databundle: false` after the first model run when all the needed data will be successfully extracted to avoid data loss.

Snakemake
===========================

Snakemake is a workflow management tool inherited by PyPSA-Earth from PyPSA-Eur. Snakemake decomposes a large software process into a set of subtasks, or ’rules’, that are automatically chained to obtain the desired output.

.. note::
  ``Snakemake``, which is one of the major dependencies, will be automatically installed in the environment pypsa-earth, thereby there is no need to install it manually.

The snakemake included in the conda environment pypsa-earth can be used to execute any custom rule with the following command:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake < your custom rule >  

Starting with essential usability features, the implemented PyPSA-Earth `Snakemake procedure <https://github.com/pypsa-meets-earth/pypsa-earth/blob/main/Snakefile>`_ that allows to flexibly execute the entire workflow with various options without writing a single line of python code. For instance, you can model the world energy system or any subset of countries only using the required data. Wildcards, which are special generic keys that can assume multiple values depending on the configuration options, help to execute large workflows with parameter sweeps and various options.


.. TODO Add Snakemake tutorial links    
