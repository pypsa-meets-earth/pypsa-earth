..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _installation:

##########################################
Installation
##########################################

The subsequently described installation steps are demonstrated as shell commands, where the path before the ``%`` sign denotes the
directory in which the commands following the ``%`` should be entered.


Clone the Repository
====================

First of all, clone the `PyPSA-Earth repository <https://github.com/pypsa-meets-earth/pypsa-earth/>`_ using the version control system ``git``.
The path to the directory into which the ``git repository`` is cloned, must **not** have any spaces!
If you do not have ``git`` installed, follow installation instructions `here <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_.

.. code:: bash

    /some/other/path % cd /some/path/without/spaces

    /some/path/without/spaces % git clone https://github.com/pypsa-meets-earth/pypsa-earth.git


.. _deps:

Install Dependencies
===============================

Python Dependencies
--------------------------------

PyPSA-Earth relies on a set of other Python packages to function.

The python package requirements are curated in the `envs/environment.yaml <https://github.com/pypsa-meets-earth/pypsa-earth/blob/main/envs/environment.yaml>`_ file.


The python package requirements are curated in the envs/environment.yaml file.
We install only `mamba` in the conda base environment to accelerate the installation.
**Please keep the base environment always clean, meaning don't install anything there!**
The environment can be installed in about 5-15 minutes and activated like this:

.. code:: bash

    .../pypsa-earth (base) % conda install -c conda-forge mamba

    .../pypsa-earth (base) % mamba env create -f envs/environment.yaml

    .../pypsa-earth (pypsa-earth) % conda activate pypsa-earth

Note please that activation is local to the currently open shell. Every time you 
open a new terminal window `pypsa-earth` environment shold be activated again to supply the workflow with all the dependencies it needs.    

In case mamba did not work for you, you might want to try conda instead:

.. code:: bash

    .../pypsa-earth (base) % conda env create -f envs/environment.yaml

    .../pypsa-earth (pypsa-earth) % conda activate pypsa-earth


For more information on conda installations you can look into :ref:`software_hints`.

Java Installation 
---------------------------------

Verify or install a java redistribution from the `official website <https://www.oracle.com/java/technologies/downloads/>`_ or equivalent. In Linux and Mac OS that is possible through the following command.

.. code:: bash

    .../pypsa-earth (pypsa-earth) % conda install -c conda-forge openjdk

To verify the successfull installation the following code can be tested from bash:

  .. code:: bash

    .../pypsa-earth % java -version

   The expected output should resemble the following:
   
   .. code:: bash
      java version "1.8.0_341"
      Java(TM) SE Runtime Environment (build 1.8.0_341-b10)
      Java HotSpot(TM) 64-Bit Server VM (build 25.341-b10, mixed mode)

Solver Installation 
---------------------------------

An optimization solver is needed to solve the mathematical problem that is build with the automated workflow.
With the goal of supporting completely open source initiative, we focus on relying on Open-Source solvers, such as 

* `CBC <https://projects.coin-or.org/Cbc>`_; 

* `GLPK <https://www.gnu.org/software/glpk/>`_;

* `WinGLPK <http://winglpk.sourceforge.net/>`_; 

* `HiGHS <https://github.com/ERGO-Code/HiGHS>`_.

To further improve performances, commercial solvers like 

* `Gurobi <http://www.gurobi.com/>`_;

* `CPLEX <https://www.ibm.com/analytics/cplex-optimizer>`_
  
(both commercial licenses with free academic options) can also be used. 

A recommended instruction to install the HiGHS solver is given `here <https://github.com/PyPSA/PyPSA/blob/633669d3f940ea256fb0a2313c7a499cbe0122a5/pypsa/linopt.py#L608-L632>`_.

Set Configuration File
================================

PyPSA-Earth has several configuration options that must be specified in a ``config.yaml`` file located in the root directory. An example configuration ``config.default.yaml`` is maintained in the repository. More details on the configuration options are in :ref:`config`.

Before first use, create a ``config.yaml`` by copying the example.

.. code:: bash

    .../pypsa-earth % cp config.default.yaml config.yaml

It makes sense to regularly check their own ``config.yaml`` against changes in the ``config.default.yaml`` when pulling a new version from the remote repository.


The snakemake included in the conda environment pypsa-earth installed with the above-mentioned procedure can be executed with the following procedure:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % .../pypsa-earth % conda activate pypsa-earth

    .../pypsa-earth (pypsa-earth) % snakemake < any command here >

