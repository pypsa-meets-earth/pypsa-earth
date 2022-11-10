..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _short_tutorial:


##########################################
Short Tutorial
##########################################

The installation procedure installs PyPSA-Earth model with all the software dependencies needed to build and run it.
To properly model any region of the Earth, PyPSA-Earth needs to download and fetch different datasets.
This section explains how to perform this data management.

How to build the tutorial model?
-------------------------------------------------

The user can explore the majority of the model's functions on a local machine by running the tutorial,
which uses fewer computational resources than the entire model does. A tutorial data kit was developed to facilitate exploring the model.
You can build it using the tutorial configuration file `config.tutorial.yaml` (placed in the project folder `pypsa-earth`).
To do that, you may want to do a reserve copy of your current configuration file and then overwrite it by a tutorial configuration:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % cp config.tutorial.yaml config.yaml


In the configuration file `config.yaml` there is a flag `retrieve_databundle` which triggers data loading and a `tutorial` flag which
determines that the loaded data belong to the light tutorial kit:

.. code:: yaml

    tutorial: true
    ...
    retrieve_databundle: true

It's recommended to set `retrieve_databundle: true` when building the model for the first time to download all needed common data files.
When the first run is completed and all the necessary data are extracted, it may be a good idea to set `retrieve_databundle: false` to avoid data loss.

How to run the model?
-------------------------------------------------

After configuration set-up, the model is ready to be built and run.
Before to actually run the workflow you may check how it will look by using `--dryrun` or `-n` Snakemake option:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 solve_all_networks --dryrun

To run the whole modeling workflow you just need the following command:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 solve_all_networks

.. TODO Explain settings of the tutorial case

This command will trigger loading of the whole dataset needed to build the model for a tutorial case if
both `tutorial` and `retrieve_databundle` flags are on. The tutorial model will run simulation of power systems in Nigeria and Benin.
Note that data load will need about 1.6Gb and model building will take a while (about 20..50 minutes).


How to analyse the solved networks?
-------------------------------------------------

The solved networks can be analysed just like any other PyPSA network (e.g. in Jupyter Notebooks).

.. code:: python

    import pypsa
    network = pypsa.Network("results/networks/elec_s_6_ec_lcopt_Co2L-4H.nc")    

We prepared an example notebook of how you could explore the above generated tutorial network in `sample-network-analysis.ipynb`.
You can find the notebook in our `notebooks repository <https://github.com/pypsa-meets-earth/documentation/tree/main/notebooks>`_.
For further inspiration how how to explore networks and also revealing the capabilities of PYPSA you may want to see the
`examples section in the PyPSA framework documentation <https://pypsa.readthedocs.io/en/latest/examples-basic.html>`_.


After playing with the tutorial model and before playing with different fucntions,
it's important to clean-up data in your model folder before to proceed further to avoid data conflicts.
You may use the `clean` rule for making so:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 clean

Generally, it's a good idea to repeat the cleaning procedure every time when the underlying data are changed.

.. note::

  This tutorial only covers Nigeria. To make the workflow run on other regions you need to use the ``config.default.yaml`` as ``config.yaml``.
  To use the model in and outside Africa, you should also read
  `How to create a model for you region of interest with PyPSA-Earth? <https://github.com/pypsa-meets-earth/pypsa-earth/discussions/505>`_
