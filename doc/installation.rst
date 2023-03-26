.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _installation:

##########################################
Installation
##########################################

The subsequently described installation steps are demonstrated as shell commands, where the path before the ``%`` sign denotes the directory in which the commands following the ``%`` should be entered.


Clone the Repository
====================
.. note::

  In order to work with the provided Jupyter notebooks in the `documentation repository <https://github.com/pypsa-meets-earth/documentation>`_, it is recommended to follow the folder structure suggested in :ref:`notebooks`.

First of all, clone the `PyPSA-Earth repository <https://github.com/pypsa-meets-earth/pypsa-earth/>`_ using the version control system ``git``.
The path to the directory into which the ``git repository`` is cloned, must **not** have any spaces.
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

The python package requirements are curated in the `envs/environment.yaml <https://github.com/pypsa-meets-earth/pypsa-earth/blob/main/envs/environment.yaml>`_ file. We install only `mamba` in the conda base environment to accelerate the installation.
**Please keep the base environment always clean, meaning don't install anything there!** It will allow to ensure compatibility of all the packages needed to work with PyPSA-Earth model.

The environment can be installed and activated like this:

.. code:: bash

    .../pypsa-earth (base) % conda install -c conda-forge mamba

    .../pypsa-earth % mamba env create -f envs/environment.yaml

    .../pypsa-earth (pypsa-earth) % conda activate pypsa-earth

Environment installation with mamba usually takes about 10-20 minutes. Note please that activation is local to the currently open shell. Every time you 
open a new terminal window, `pypsa-earth` environment should be activated again to supply the workflow with all the dependencies it needs.    

In case mamba did not work for you, you might want to try conda instead:

.. code:: bash

    .../pypsa-earth % conda env create -f envs/environment.yaml

    .../pypsa-earth (pypsa-earth) % conda activate pypsa-earth


For more on information on how to install conda and work with it you can look into :ref:`software_hints`.

Java Installation 
---------------------------------

PyPSA-Earth currently needs Java redistribution to work properly. To check if Java is still installed you can request it's version from a terminal:

  .. code:: bash

    .../pypsa-earth % java --version

The expected output should resemble the following:
   
   .. code:: bash
      java version "1.8.0_341"
      Java(TM) SE Runtime Environment (build 1.8.0_341-b10)
      Java HotSpot(TM) 64-Bit Server VM (build 25.341-b10, mixed mode)

In case you don't have Java, you have to install it from the `official website <https://www.oracle.com/java/technologies/downloads/>`_ or equivalent. In Linux and Mac OS that is possible through the following command:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % conda install -c conda-forge openjdk


Solver Installation 
---------------------------------

An optimization solver is needed to solve the mathematical problem that is build with the automated workflow of PyPSA-Earth.
With the goal of supporting completely open source initiative, we focus on relying on Open-Source solvers, such as 

* `CBC <https://projects.coin-or.org/Cbc>`_; 

* `GLPK <https://www.gnu.org/software/glpk/>`_ and `WinGLPK <http://winglpk.sourceforge.net/>`_ (is included into pypsa-earth environment and installed automatically during environment creation); 

* `HiGHS <https://github.com/ERGO-Code/HiGHS>`_.

To further improve performances, commercial solvers like 

* `Gurobi <http://www.gurobi.com/>`_;

* `CPLEX <https://www.ibm.com/analytics/cplex-optimizer>`_.
  
(both commercial licenses with free academic options) can also be used. 

A recommended instruction to install the HiGHS solver is given `here <https://github.com/PyPSA/PyPSA/blob/633669d3f940ea256fb0a2313c7a499cbe0122a5/pypsa/linopt.py#L608-L632>`_.

Set Configuration File
================================

PyPSA-Earth has several configuration options that must be specified in a ``config.yaml`` file located in the project directory. An example configuration ``config.default.yaml`` is maintained in the repository. More details on the configuration options are in :ref:`config` section.

Before first use, create a ``config.yaml`` by copying the example.

.. code:: bash

    .../pypsa-earth % cp config.default.yaml config.yaml

It makes sense to regularly check their own ``config.yaml`` against changes in the ``config.default.yaml`` when pulling a new version from the remote repository.

Install Jupyter Lab
================================

We use Jupyter notebooks to share examples on how to use the model and analyse the results. VSCode supports working with Jupyter Notebooks natively. In case you are using different IDE and don't have Jupyter notebooks pre-installed you can install jupyter lab (new jupyter notebooks) with the `ipython kernel installation <http://echrislynch.com/2019/02/01/adding-an-environment-to-jupyter-notebooks/>`_ and test if your jupyter lab works:

.. code:: bash

    .../pypsa-earth % ipython kernel install --user --name=pypsa-earth
    .../pypsa-earth % jupyter lab
