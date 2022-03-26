..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Africa authors

  SPDX-License-Identifier: CC-BY-4.0

.. _introduction:

##########################################
Introduction
##########################################

PyPSA meets Africa aims at providing an open-source energy system model, named PyPSA-Africa, that describes the African energy system **by the end of 2021**.

.. raw:: html

    <iframe width="832" height="468" src="https://www.youtube.com/embed/E0V0T4U9nmQ" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>


The Mission
===========

The PyPSA Africa team foreseas the use of Open Energy Modelling for Africa as a disruptive approach to support policy makers, utilities,
developers and researchers in performing the most affordable, reliable, transparent, reproducible and informed decisions that Africa deserves.


We are committed to give tools that can efficiently provide instruments to all stakeholders to perform the best informed decisions and promote coordinated
planning and dispatch to maximize the efficient use of energy that can sustain stable sustainable growth.

..
    Despite being home to almost twice the population of Europe, Africa energy demand is a quarter of Europe [OWD]_.
    Access to energy is also very diverse and, according to IEA, still around 600mln people have no access to electricity in Sub-Saharan Africa [WEO2021]_.

    .. [OWD] https://ourworldindata.org/grapher/primary-energy-consumption-by-region
    .. [WEO2021] https://www.iea.org/data-and-statistics/charts/people-without-access-to-electricity-in-sub-saharan-africa-2000-2021

Project structure
=================

In PyPSA Africa, six main work packages are in place to tackle the continent-wide challenges:

- **WP1. Demand modelling**: to characterize the demand profiles for each location both in space and time resolution
- **WP2. Conventional generator modelling**: define the type of generators that are installed in each substation
- **WP3. RES modelling**: to characterize the available renewable production in space and time resolution, given technical, social, and land constraints
- **WP4. Land coverage constraint modelling**: to identify limits and constraints where assets can or cannot be installed
- **WP5. Network and substation modelling**: drawing the model of the network to feed the PyPSA modeling suite
- **WP6. Data creation and validation**: when data are scarce,
  by using high resolution maps we aim at filling the gaps and update the datasets. The AI detection is tackled in the package
  `detect_energy <https://github.com/pypsa-meets-africa/detect_energy>`_


For more details, read our preliminary conference paper (more to come):

- D. Kirli et al., "PyPSA meets Africa: Developing an open source electricity network model of the African continent," 2021 IEEE AFRICON, 2021, pp. 1-6, doi: `10.1109/AFRICON51333.2021.9570911 <https://doi.org/10.1109/AFRICON51333.2021.9570911>`_.

Please use the following BibTeX: ::

    @INPROCEEDINGS{9570911,
        author={Kirli, Desen and Hampp, Johannes and van Greevenbroek, Koen and Grant, Rebecca and Mahmood, Matin and Parzen, Maximilian and Kiprakis, Aristides},
        booktitle={2021 IEEE AFRICON},
        title={PyPSA meets Africa: Developing an open source electricity network model of the African continent},
        year={2021},
        volume={},
        number={},
        pages={1-6},
        doi={10.1109/AFRICON51333.2021.9570911}
    }


Workflow
========

The data workflow is entirely managed through the `Snakemake <https://snakemake.bitbucket.io/>`_ workflow management system that ease the development of automated reproducible workflow procedure.
Based on the configuration files, ``snakefile`` iterate the execution of the scripts in the ``scripts`` folder, according to pre-defined `rules`, which guide what are the inputs and outputs for each script.
Then, once ``Snakefile`` tool is executed, it automatically identify how to order the executions of the rules, according to the pipeline of dependencies of inputs and outputs.

For example, by executing the following snakemake routine

.. code:: bash

    .../pypsa-africa % snakemake -j 1 networks/elec_s_128.nc

the following workflow is automatically executed.

.. image:: img/workflow_introduction.png
    :align: center

Each block represents the execution of a script in the ``script`` directory and the arrows describe the input/output dependencies between the rules.
Note that none of the above has been manually performed, but ``snakemke`` automatically recognizes the correct order and executes the rules accordingly.

.. note::
    For reproducibility purposes, the image can be obtained through
    ``snakemake --dag networks/elec_s_128.nc | dot -Tpng -o workflow.png``
    using `Graphviz <https://graphviz.org/>`_



Folder structure
================

The content in this package is organized in folders as described below; for more details, please see the documentation.

- ``data``: Includes input data that is not produced by any ``snakemake`` rule.
- ``scripts``: Includes all the Python scripts executed by the ``snakemake`` rules.
- ``notebooks``: Stores useful notebooks to investigate the raw data, the intermediate data and the results, as well as to generate plots.
- ``resources``: Stores intermediate results of the workflow which can be picked up again by subsequent rules.
- ``networks``: Stores intermediate, unsolved stages of the PyPSA network that describes the energy system model.
- ``results``: Stores the solved PyPSA network data, summary files and plots.
- ``benchmarks``: Stores ``snakemake`` benchmarks.
- ``logs``: Stores log files about solving, including the solver output, console output and the output of a memory logger.
- ``envs``: Stores the conda environment files to successfully run the workflow.


Licence
=======

PyPSA-Africa work is released under multiple licenses:

* All original source code is licensed as free software under `GPL-3.0 License <https://github.com/pypsa-meets-africa/pypsa-africa/blob/main/LICENSE>`_.
* The documentation is licensed under `CC-BY-4.0 <https://creativecommons.org/licenses/by/4.0/>`_.
* Configuration files are mostly licensed under `CC0-1.0 <https://creativecommons.org/publicdomain/zero/1.0/>`_.
* Data files are licensed under different licenses as noted below.

Licenses and urls of the data used in PyPSA-Africa:
.. csv-table::
   :header-rows: 1
   :file: configtables/licenses.csv

* *BY: Attribute Source*
* *NC: Non-Commercial Use Only*
* *SA: Share Alike*
