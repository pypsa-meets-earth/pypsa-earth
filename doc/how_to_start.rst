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

`clean` to avoid data conflicts (it's a good idea to repeat it when the underlying data are changed)

`retrieve_databundle: true`

Execute the following code on the shell to download initial files. Please, note that around **20Gb zipped files will be downloaded**, make sure you have a stable connection, time and around 50 Gb available in your system. If no errors show up, then you can proceed.

.. code:: bash
    snakemake -j 1 `solve_all_networks`

`retrieve_databundle: false` to avoid data loss 

Snakemake
===========================

- Snakemake

.. note::
  ``Snakemake``, which is one of the major dependencies, will be automatically installed in the environment pypsa-earth, thereby there is no need to install it manually.

The snakemake included in the conda environment pypsa-earth installed with the above-mentioned procedure can be executed with the following procedure:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake < any command here >  

.. TODO Add Snakemake tutorial links    