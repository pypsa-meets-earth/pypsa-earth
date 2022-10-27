..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _how_to_start:


##########################################
How to run the model
##########################################

Data are needed to run a model.

Data workflow management is one of the core PyPSA-Earth features.

Tutorial run
------------------------------

The entire distribution, including the data for most parts on Earth, is very heavy (>40Gb for Africa) and it involves a large number of files. To simplify the installation of the github folder, the main source code is available in the Github folder, whereas the data are stored in cloud.

`tutorial` flag in the `config.yaml`

.. code:: bash
    snakemake -j 1 `solve_all_networks`

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