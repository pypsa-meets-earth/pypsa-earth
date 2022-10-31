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

A tutorial data kit was developed to facilitate exploring the model. You can use it by using the tutorial configuration file `config.tutorial.yaml` (placed in the project folder `pypsa-earth`). To do that, you may want to do a reserve copy of your current configuration file and then overwrite it by a tutorial configuration:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % cp config.tutorial.yaml config.yaml


In the configuration file `config.yaml` there is a flag `retrieve_databundle` which triggers data loading and a `tutorial` flag which determines that the loaded data belong to the light tutorial kit:

.. code:: yaml

    tutorial: true
    ...
    retrieve_databundle: true

It's recommended to set `retrieve_databundle: true` when building the model first time to download all the common data files needed. When the first run is completed and all the necessary data extracted, it may be a good idea to set `retrieve_databundle: false` to avoid data loss.

How to run the model?
------------------------------

After configuration set-up, you just need run the modeling workflow with the following command:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 solve_all_networks

.. TODO Explain settings of the tutorial case

This command will trigger loading of the whole dataset needed to build the model for a tutorial case if both `tutorial` and `retrieve_databundle` flags are on. The tutorial model will run simulation of power systems in Nigeria and Benin. Note please that data load will need about 1.6Gb and model building will take a while (about 20..50 minutes).

Before to actually run the workflow you may check how it will look by using `dryrun` Snakemake option:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 solve_all_networks --dryrun


How to analyse the solved networks?
------------------------------

The solved networks can be analysed just like any other PyPSA network (e.g. in Jupyter Notebooks).

.. code:: python

    import pypsa
    network = pypsa.Network("results/networks/elec_s_6_ec_lcopt_Co2L-4H.nc")    

For inspiration, you may want to have a look on the `examples section in the PyPSA documentation <https://pypsa.readthedocs.io/en/latest/examples-basic.html>`_.

After playing with the tutorial model, it's important to clean-up data in your model folder before to proceed further to avoid data conflicts. You may use the `clean` rule for making so:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 clean

Generally, it's a good idea to repeat the cleaning procedure every time when the underlying data are changed.

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
