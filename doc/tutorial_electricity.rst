.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _tutorial_electricity:


##########################################
Tutorial: Electricity
##########################################

.. note::
    If you have not yet installed PyPSA-Earth, please refer to the :ref:`installation` section.

To properly model any region of the Earth, it is first crucial to get familiar with a tutorial where a simpler model is considered.
This section explains how to run and analyze the tutorial model.

Build the tutorial electricity-only model
---------------------------

The user can explore the majority of the model's functions on a local machine by running the tutorial, which uses fewer computational resources than the entire model does. A tutorial data kit was developed to facilitate exploring the model.
You can build it using the tutorial configuration file ``config.tutorial.yaml`` (placed in the project folder ``pypsa-earth``). It includes parts deviating from the default config file ``config.default.yaml``, which are necessary to run the tutorial model.
By default, PyPSA-Earth reads configuration parameters of simulation from ``config.yaml`` file located in ``pypsa-earth`` folder.
Thus, to run the tutorial model, ``config.tutorial.yaml`` needs to be stored as ``config.yaml``:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % cp config.tutorial.yaml config.yaml

.. note::
    You may want to do a reserve copy of your current configuration file (``config.yaml``) as it will be overwritten by a tutorial configuration.

In the configuration file ``config.yaml`` there is a flag ``retrieve_databundle`` which triggers data loading and a ``tutorial`` flag which determines that the loaded data belong to the tutorial kit. Currently the tutorial can be run only for Nigeria ("NG"), Benin ("BJ"), Botswana ("BW") and Morocco ("MA").

.. code:: yaml

    tutorial: true
    ...
    enable:
        retrieve_databundle: true

.. note::
    It's recommended to set ``retrieve_databundle: true`` when building the model for the first time to download all needed common data files.
    When the first run is completed and all the necessary data are extracted, it may be a good idea to set ``retrieve_databundle: false`` to avoid data loss.

Run the model
---------------------

After configuration set-up, the model is ready to be built and run.
Before running the workflow you may check how it will look by using `--dryrun` or `-n` Snakemake option:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 solve_all_networks --dryrun

To run the whole model workflow you just need the following command:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 solve_all_networks

.. TODO Explain settings of the tutorial case

This command will trigger loading of the whole dataset needed to build the model for a tutorial case if
both ``tutorial`` and ``retrieve_databundle`` flags are on. The tutorial model will run simulation of power systems in Nigeria and Benin.
Note that data load will need about 1.6GB and model building will take a while (about 20-50 minutes).


Analyse the solved networks
------------------------------------

The solved networks can be analysed just like any other PyPSA network (e.g. in Jupyter Notebooks).

.. code:: python

    import pypsa

    network = pypsa.Network("results/networks/elec_s_6_ec_lcopt_Co2L-4H.nc")

The video below shows how to analyse solved PyPSA-Eur networks in Jupyter Notebooks.
Fabian Neumann did a great job explaining the basics of PyPSA and how to use it for analysis.

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/mAwhQnNRIvs" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

We also prepared an example notebook such that you can explore the tutorial network yourself.
Just open in our `notebooks repository <https://github.com/pypsa-meets-earth/documentation/tree/main/notebooks>`_
the file ``sample-network-analysis.ipynb``. For further inspiration on what you can analyse and do with PyPSA,
you can explore the `examples section in the PyPSA framework documentation <https://pypsa.readthedocs.io/en/latest/examples-basic.html>`_.

After playing with the tutorial model and before playing with different functions,
it's important to clean-up data in your model folder before to proceed further to avoid data conflicts.
You may use the ``clean`` rule for making so:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 clean

Generally, it's a good idea to repeat the cleaning procedure every time when the underlying data are changed to avoid conflicts between run settings corresponding to different scenarios.

It is also possible to make manual clean-up removing folders "resources", "networks" and "results". Those folders store the intermediate output of the workflow and if you don't need them anymore it is safe to delete them.

.. note::

  This tutorial only covers Nigeria and Benin. To make the workflow run on other regions you need to use the ``config.default.yaml`` as ``config.yaml``.
  To use the model in and outside Africa, you should also read
  `How to create a model for you region of interest with PyPSA-Earth? <https://github.com/pypsa-meets-earth/pypsa-earth/discussions/505>`_

:ref:`tutorial` section elaborates on building and running a full PyPSA-Earth model.
