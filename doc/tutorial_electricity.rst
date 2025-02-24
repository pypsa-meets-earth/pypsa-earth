.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _tutorial_electricity:


##########################################
Tutorial: Electricity
##########################################

.. note::
    If you have not yet installed PyPSA-Earth, please refer to the :ref:`installation` section.

To properly model any region of the Earth, it is first crucial to get familiar with a tutorial
where a simpler model is considered. This section explains how to run and analyze the tutorial model.

Build the tutorial electricity-only model
---------------------------------------------

The user can explore the majority of the model's functions on a local machine by running the tutorial,
which uses fewer computational resources than the entire model does. A tutorial data kit was developed
to facilitate exploring the model.
You can build it using the tutorial configuration file ``config.tutorial.yaml`` (placed in the project
folder ``pypsa-earth``). It includes parts deviating from the default config file ``config.default.yaml``,
which are necessary to run the tutorial model. By default, PyPSA-Earth reads configuration parameters
of simulation from ``config.yaml`` file located in ``pypsa-earth`` folder. Thus, to run the tutorial
model, ``config.tutorial.yaml`` needs to be stored as ``config.yaml``:


How to configure runs for the tutorial model
---------------------------------------------

The model can be adapted to include any selected country. But this tutorial is limited to ``Nigeria ("NG")``,
``Benin ("BJ")``, ``Botswana ("BW")`` and ``Morocco ("MA")``.

.. code:: yaml

    countries: ["NG", "BJ"]

It's recommended to set ``retrieve_databundle: true`` when building the model for the first time to download all needed common data files.
When the first run is completed and all the necessary data are extracted, it may be a good idea to set ``retrieve_databundle: false`` to avoid data loss.

.. code:: yaml

    enable:
        retrieve_databundle: true

The scenario is defined by the number of clusters and the optimization options. The tutorial model
is set to have 6 clusters and the optimization option "Co2L-4H" which translates absolute carbon-dioxide
emission limit with the model resampled to 4H resolution.

.. code:: yaml

    scenario:
        clusters: [6]
        opts: [Co2L-4H]

The temporal scope is set to a single week. This is to make sure that the model completes in no time.

.. code:: yaml

    snapshots:
        start: "2013-03-1"
        end: "2013-03-7"

.. note::
    For more information on the configuration file, please refer to the :ref:`config` section.

Run the model
---------------------

After configuration set-up, the model is ready to be built and run.
Open a terminal, go into the PyPSA-Earth directory, and activate the pypsa-earth environment with

.. code:: bash

    .../pypsa-earth $ conda activate pypsa-earth

You then need to copy the tutorial config file to ``config.yaml``

.. code:: bash

    .../pypsa-earth (pypsa-earth) $ cp config.tutorial.yaml config.yaml

.. note::
    If you previously have a ``config.yaml file``, You may want to reserve a copy of
    your current configuration file (``config.yaml``) as it will be overwritten by a tutorial configuration.

Before running the workflow you may check how it will look by using ``--dryrun`` or ``-n`` Snakemake option:

.. code:: bash

    .../pypsa-earth (pypsa-earth) $ snakemake -j 1 solve_all_networks --dryrun

This triggers a workflow of multiple preceding jobs that depend on each rule's inputs and outputs:

.. graphviz::
    :class: full-width
    :align: center

        digraph snakemake_dag {
        graph[bgcolor=white, margin=0];
        node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
        edge[penwidth=2, color=grey];
        0[label = "solve_all_networks", color = "0.33 0.6 0.85", style="rounded"];
        1[label = "solve_network", color = "0.37 0.6 0.85", style="rounded"];
        2[label = "prepare_network\nll: copt\nopts: Co2L-4H", color = "0.02 0.6 0.85", style="rounded"];
        3[label = "add_extra_components", color = "0.50 0.6 0.85", style="rounded"];
        4[label = "cluster_network\nclusters: 6", color = "0.31 0.6 0.85", style="rounded"];
        5[label = "simplify_network\nsimpl: ", color = "0.23 0.6 0.85", style="rounded"];
        6[label = "add_electricity", color = "0.43 0.6 0.85", style="rounded"];
        7[label = "build_renewable_profiles\ntechnology: onwind", color = "0.20 0.6 0.85", style="rounded"];
        8[label = "build_natura_raster", color = "0.07 0.6 0.85", style="rounded"];
        9[label = "retrieve_databundle_light", color = "0.04 0.6 0.85", style="rounded"];
        10[label = "build_shapes", color = "0.59 0.6 0.85", style="rounded"];
        11[label = "build_powerplants", color = "0.38 0.6 0.85", style="rounded"];
        12[label = "base_network", color = "0.47 0.6 0.85", style="rounded"];
        13[label = "build_osm_network", color = "0.52 0.6 0.85", style="rounded"];
        14[label = "clean_osm_data", color = "0.13 0.6 0.85", style="rounded"];
        15[label = "download_osm_data", color = "0.32 0.6 0.85", style="rounded"];
        16[label = "build_bus_regions", color = "0.06 0.6 0.85", style="rounded"];
        17[label = "build_renewable_profiles\ntechnology: offwind-ac", color = "0.20 0.6 0.85", style="rounded"];
        18[label = "build_renewable_profiles\ntechnology: offwind-dc", color = "0.20 0.6 0.85", style="rounded"];
        19[label = "build_renewable_profiles\ntechnology: solar", color = "0.20 0.6 0.85", style="rounded"];
        20[label = "build_renewable_profiles\ntechnology: hydro", color = "0.20 0.6 0.85", style="rounded"];
        21[label = "retrieve_cost_data\nyear: 2030", color = "0.44 0.6 0.85", style="rounded"];
        22[label = "build_demand_profiles", color = "0.51 0.6 0.85", style="rounded"];
        1 -> 0
        2 -> 1
        3 -> 2
        21 -> 2
        4 -> 3
        21 -> 3
        5 -> 4
        10 -> 4
        21 -> 4
        6 -> 5
        21 -> 5
        16 -> 5
        10 -> 5
        7 -> 6
        17 -> 6
        18 -> 6
        19 -> 6
        20 -> 6
        12 -> 6
        21 -> 6
        11 -> 6
        10 -> 6
        22 -> 6
        8 -> 7
        9 -> 7
        10 -> 7
        11 -> 7
        16 -> 7
        9 -> 8
        9 -> 10
        12 -> 11
        14 -> 11
        10 -> 11
        13 -> 12
        10 -> 12
        14 -> 13
        10 -> 13
        15 -> 14
        10 -> 14
        10 -> 16
        12 -> 16
        8 -> 17
        9 -> 17
        10 -> 17
        11 -> 17
        16 -> 17
        8 -> 18
        9 -> 18
        10 -> 18
        11 -> 18
        16 -> 18
        8 -> 19
        9 -> 19
        10 -> 19
        11 -> 19
        16 -> 19
        8 -> 20
        9 -> 20
        10 -> 20
        11 -> 20
        16 -> 20
        12 -> 22
        16 -> 22
        9 -> 22
        10 -> 22
    }

In the terminal, this will show up as a list of jobs to be run:

.. code:: console

    Building DAG of jobs...
    Job stats:
    job                          count
    -------------------------  -------
    add_electricity                  1
    add_extra_components             1
    base_network                     1
    build_bus_regions                1
    build_demand_profiles            1
    build_natura_raster              1
    build_osm_network                1
    build_powerplants                1
    build_renewable_profiles         5
    build_shapes                     1
    clean_osm_data                   1
    cluster_network                  1
    download_osm_data                1
    prepare_network                  1
    retrieve_cost_data               1
    retrieve_databundle_light        1
    simplify_network                 1
    solve_all_networks               1
    solve_network                    1
    total                           23


To run the whole model workflow you just need the following command:

.. code:: bash

    .../pypsa-earth (pypsa-earth) $  snakemake -j 1 solve_all_networks

You can also run the tutorial model using the tutorial config directly by using the following command:

.. code:: bash

    .../pypsa-earth (pypsa-earth) $ snakemake -j 1 solve_all_networks --configfile config.tutorial.yaml


This command will trigger loading of the whole dataset needed to build the model for a tutorial case if
both ``tutorial`` and ``retrieve_databundle`` flags are on. The tutorial model will run simulation of power systems in Nigeria and Benin.
Note that data load will need about 1.6GB and model building will take a while (about 20-50 minutes).

.. note::

    It is good practice to perform a dry-run using the option -n, before you commit to a run:

    .. code:: bash

        .../pypsa-earth (pypsa-earth) $ snakemake solve_all_networks -n



Analyse the solved networks
------------------------------------

The solved networks can be analysed just like any other PyPSA network (e.g. in Jupyter Notebooks).

.. code:: python

    import pypsa

    network = pypsa.Network("results/networks/elec_s_6_ec_lcopt_Co2L-4H.nc")

The video below shows how to analyse solved PyPSA-Eur networks in Jupyter Notebooks.
Fabian Neumann did a great job explaining the basics of PyPSA and how to use it for analysis.

.. raw:: html

    <iframe width="832" height="468" src="https://www.youtube.com/embed/mAwhQnNRIvs" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

We also prepared an example notebook such that you can explore the tutorial network yourself.
Just open in our `notebooks repository <https://github.com/pypsa-meets-earth/documentation/tree/main/notebooks>`_
the file ``sample-network-analysis.ipynb``. For further inspiration on what you can analyse and do with PyPSA,
you can explore the `examples section in the PyPSA framework documentation <https://pypsa.readthedocs.io/en/latest/getting-started/quick-start.html>`_.

After playing with the tutorial model and before playing with different functions,
it's important to clean-up data in your model folder before to proceed further to avoid data conflicts.
You may use the ``clean`` rule for making so:

.. code:: bash

    .../pypsa-earth (pypsa-earth) $ snakemake -j 1 clean

Generally, it's a good idea to repeat the cleaning procedure every time when the underlying data are changed to avoid conflicts between run settings corresponding to different scenarios.

It is also possible to make manual clean-up removing folders "resources", "networks" and "results". Those folders store the intermediate output of the workflow and if you don't need them anymore it is safe to delete them.

.. note::

  This tutorial only covers Nigeria and Benin. To make the workflow run on other regions you need to use the ``config.default.yaml`` as ``config.yaml``.
  To use the model in and outside Africa, you should also read
  `How to create a model for you region of interest with PyPSA-Earth? <https://github.com/pypsa-meets-earth/pypsa-earth/discussions/505>`_

:ref:`tutorial` section elaborates on building and running a full PyPSA-Earth model.
