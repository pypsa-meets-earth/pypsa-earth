.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _customization_run:

####################
2. Running the model
####################

To solve full optimization problem, it is important to pick a solver in `config.yaml` file. For instance, this tutorial uses the open-source solver glpk and does not rely
on the commercial solvers such as Gurobi or CPLEX (for which free academic licenses are available).

.. code:: yaml

    solver:
        name: glpk

.. note::

    ``glpk`` can solve the network with low temporal and spacial resolution. To make a full model run, it is adviced to use ``CPLEX``, ``Gurobi``, or open-source `HIGHs <https://highs.dev/>`__.

To execute a full model run, the following command needs to be applied:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake --cores 1 solve_all_networks

Here, ``snakemake`` is a workflow management tool inherited by PyPSA-Earth from PyPSA-Eur.
Snakemake decomposes a large software process into a set of subtasks, or ’rules’, that are automatically chained to obtain the desired output.
The flag ``--cores 1`` dictates the number of CPU cores allocated for the process. Notably, the ``solve_all_network`` rule within Snakemake orchestrates the process of solving the network. 

.. note::

  ``Snakemake``, which is one of the major dependencies, will be automatically installed in ``pypsa-earth`` environment, thereby there is no need to install it manually.

The snakemake included in the ``pypsa-earth`` conda environment pypsa-earth can be used to execute any custom rule with the following command:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake < your custom rule >

Starting with essential usability features, PyPSA-Earth `Snakemake procedure <https://github.com/pypsa-meets-earth/pypsa-earth/blob/main/Snakefile>`_ enables users to execute the entire workflow flexibly with diverse options, requiring no Python coding. For example, users can model the global energy system or specific subsets of countries using only necessary data. Wildcards, serving as special generic keys, adapt to multiple values based on configuration options, facilitating the execution of extensive workflows with parameter sweeps and diverse options.

You can execute some parts of the workflow in case you are interested in some specific parts.
E.g. power grid topology may be extracted and cleaned with the following command which refers to the script name:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 clean_osm_data

Solar profile for the requested area may be calculated using the output name:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 resources/renewable_profiles/profile_solar.nc