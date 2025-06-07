.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _installation:

##########################################
Installation
##########################################

Before installing PyPSA-Earth on your computer, it's crucial to ensure that your system meets the necessary hardware and software requirements. The following sections outline the prerequisites in terms of hardware and software. Additionally, detailed installation guidelines for required software tools will be provided, followed by step-by-step instructions for installing PyPSA-Earth.

Hardware Requirements
=====================
Ensure that your system meets the minimum hardware specifications to run PyPSA-Earth effectively. Recommended hardware specifications may include:

* 8-16 GB RAM and adequate CPU (at-least 2-cores)

* Storage (HDD/SSD) capacity depends on the region of interest. Africa model requires 40 GB, the world `--`` 250 GB, a single country `--` between 1-10 GB. Tutorial requires just below 2 GB. Thus, considering all required software tools, at least 40 GB of storage space is recommended.

.. note::

    The subsequently described installation steps are demonstrated as shell commands, where the path before the ``%`` sign denotes the directory in which the commands following the ``%`` should be entered.

Software Prerequisites
======================
Prior to installing PyPSA-Earth, you'll need to ensure the following software tools are installed and configured on your system:

* `Miniconda <https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html>`_
* `Git <https://git-scm.com/downloads>`_
* `VSCode <https://code.visualstudio.com/>`_ (or any other IDE)
* `Java <https://www.oracle.com/java/technologies/downloads/>`_

Miniconda
---------
To use packages in python, it is highly recommended to use a ``conda`` package manager, such as `miniconda <https://docs.conda.io/projects/miniconda/en/latest/>`__. You may check if ``conda`` is already installed on your system with

.. code:: bash

    $ conda --version

If ``conda`` is not installed, follow `miniconda installation guide <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_.
For more on information on how to install conda and work with it you can look into :ref:`software_hints`.

Git
---
`Git <https://git-scm.com/>`__ is a free open-source tool that facilitates tracking changes in the code development and enable to coordinate the parallel software development between many developers.
Download and install ``git`` to your system using the following `link <https://git-scm.com/downloads>`__.
It is highly recommended to `learn the git basics <https://git-scm.com/doc>`__.

VSCode
------
In order to write and debug python code, you need an Integrated Development Environment (IDE) that is a software used to write code. We recommend `Visual Studio Code <https://code.visualstudio.com/>`_, which is freely available online and provides an easy-to-use interface with Git. Obviously, any alternatives like `PyCharm <https://www.jetbrains.com/pycharm/>`_ or `Sublime <https://www.sublimetext.com/>`_ will work as well.

Java
----
PyPSA-Earth currently needs Java redistribution to work properly. To check if Java is still installed you can request it's version from a terminal:

.. code:: bash

    $ java --version

The expected output should resemble the following:

.. code:: bash

    java version "1.8.0_341"
    Java(TM) SE Runtime Environment (build 1.8.0_341-b10)
    Java HotSpot(TM) 64-Bit Server VM (build 25.341-b10, mixed mode)

In case you don't have Java, you have to install it from the `official website <https://www.oracle.com/java/technologies/downloads/>`_ or equivalent. Depending on the version of OS, download and install ``JDK 21`` or ``JDK 17`` from the link.


Installation with Conda/Mamba
========================

Clone the Repository
--------------------
.. note::

  In order to work with the provided Jupyter notebooks in the `documentation repository <https://github.com/pypsa-meets-earth/documentation>`__, it is recommended to follow the folder structure suggested in :ref:`notebooks`.

First of all, clone the `PyPSA-Earth repository <https://github.com/pypsa-meets-earth/pypsa-earth/>`__ using the version control system ``git``.
The path to the directory into which the ``git repository`` is cloned, must **not** have any spaces.
The following commands can be executed in command prompt of ``miniconda``, terminal of ``VSCode``, or in ``Git Bash``.

.. code:: bash

    $ git clone https://github.com/pypsa-meets-earth/pypsa-earth.git

Install Dependencies
-------------------------
PyPSA-Earth relies on a set of other Python packages to function.

The python package requirements are located in the `envs/environment.yaml <https://github.com/pypsa-meets-earth/pypsa-earth/blob/main/envs/environment.yaml>`_ file. We install only `mamba` in the conda base environment to accelerate the installation.
**Please keep the base environment always clean, meaning don't install anything there!** It will allow to ensure compatibility of all the packages needed to work with PyPSA-Earth model.

There are also regularly updated locked environment files for
each platform generated with conda-lock to ensure reproducibility. Choose the correct file for your platform:

* For Intel/AMD processors:

  - Linux: ``envs/linux-64.lock.yaml``

  - macOS: ``envs/osx-64.lock.yaml``

  - Windows: ``envs/win-64.lock.yaml``

* For ARM processors:

  - macOS (Apple Silicon): ``envs/osx-arm64.lock.yaml``

  - Linux (ARM): Currently not supported via lock files; requires building certain packages, such as ``PySCIPOpt``, from source

We recommend using these locked files for a stable environment.

.. note::

    you can check and verify your platform with ``conda info``

.. code:: console

    $ conda install -c conda-forge mamba

    $ mamba env create -f envs/linux-64.lock.yaml # select the appropriate file for your platform

    $ conda activate pypsa-earth

Environment installation with mamba usually takes about 10-20 minutes. Note please that activation is local to the currently open shell. Every time you
open a new terminal window, `pypsa-earth` environment should be activated again to supply the workflow with all the dependencies it needs.

In case mamba did not work for you, you might want to try conda instead:

.. code:: bash

    $ conda env create -f envs/linux-64.lock.yaml

    $ conda activate pypsa-earth


Generating the Lock Files
------------------------
If a pre-generated lock file is not available for your platform (e.g., ``aarch64``, ARM Mac, etc.), you can generate one using the following:

1. Ensure ``conda-lock`` is installed:

   .. code-block:: bash

      $ conda install conda-lock -c conda-forge

2. Generate lock files for target platforms:

   .. code-block:: bash

      $ conda-lock lock -p <your-platform> -k env -f envs/environment.yaml

For platform codes, refer to the `conda-lock documentation <https://conda.github.io/conda-lock/>`_ or use ``conda info`` to determine your platform.


For more on information on how to install conda and work with it you can look into :ref:`software_hints`.

To confirm the installation, run the following command in the activated environment:

.. code:: bash

    $ snakemake --version


Solver Installation
--------------------
An optimization solver is needed to solve the mathematical problem that is build with the automated workflow of PyPSA-Earth.
With the goal of supporting completely open source initiative, we focus on relying on Open-Source solvers, such as

* `CBC <https://projects.coin-or.org/Cbc>`_;

* `GLPK <https://www.gnu.org/software/glpk/>`_ and `WinGLPK <http://winglpk.sourceforge.net/>`_ (is included into pypsa-earth environment and installed automatically during environment creation);

* `HiGHS <https://github.com/ERGO-Code/HiGHS>`_.

To further improve performances, commercial solvers like

* `Gurobi <http://www.gurobi.com/>`_;

* `CPLEX <https://www.ibm.com/analytics/cplex-optimizer>`_.

(both commercial licenses with free academic options) can also be used.

.. note::

    No need to install ``glpk`` separately, as they are included in ``envs/environment.yaml`` and installed during ``conda`` environment creation.
    However, solving capabilities of ``glpk`` are limited.
    To run the model with high temporal and spatial resolution, it is recommended to use ``cplex``, ``gurobi``, or ``highs``.

A recommended instruction to install the HiGHS solver is given `here <https://github.com/PyPSA/PyPSA/blob/633669d3f940ea256fb0a2313c7a499cbe0122a5/pypsa/linopt.py#L608-L632>`_.


Install Jupyter Lab
================================

We use Jupyter notebooks to share examples on how to use the model and analyse the results. ``VSCode`` supports working with Jupyter Notebooks natively. In case you are using different IDE and don't have Jupyter notebooks pre-installed you can install jupyter lab (new jupyter notebooks) with the `ipython kernel installation <http://echrislynch.com/2019/02/01/adding-an-environment-to-jupyter-notebooks/>`_ and test if your jupyter lab works:

.. code:: bash

    $ ipython kernel install --user --name=pypsa-earth
    $ jupyter lab


Alternate installation with Docker
================================
This is an alternative way to create a development environment for PyPSA-Earth. This method is useful for users who are not familiar with programming or Python, or who do not want to install Python on their local machine. It uses Docker containers to create a development environment for PyPSA-Earth.

This section provides a step-by-step guide on how to set up and use Docker containers to run PyPSA-Earth.

Steps:

1. Install Docker: Follow the instructions for your operating system:

   * `Windows <https://docs.docker.com/desktop/install/windows-install/>`_
   * `Linux <https://docs.docker.com/desktop/install/linux/>`_
   * `MacOS <https://docs.docker.com/desktop/install/mac-install/>`_

.. note::

    Ensure Docker is installed on your system.

2. You can use the link `here <https://code.visualstudio.com/download>`_ to install Visual Studio Code on your operating system. Ensure to select the most compatible file for your operating system.

3. Install GitHub Desktop for your OS `here <https://desktop.github.com/download/>`_.

4. Clone the repository:

   * Open GitHub Desktop.
   * Click on "File" in the top left corner.
   * Click on "Clone Repository".
   * Paste the following URL in the URL field:

   .. code:: bash

    https://github.com/pypsa-meets-earth/pypsa-earth.git

   * Click on "Clone".
   * Choose the location where you want to save the repository.
   * Click on "Current Branch: main" and select `devContainers`.
   * Click on "Open in Visual Studio Code".

   ``The repository will be cloned to your local machine.``

5. Rebuild and open in a container:

   * Open the repository in VSCode.
   * Click on the icon in the far bottom left corner of the VSCode window.
   * Click on "Reopen in Container".
   * Wait for the container to build and open the repository in the container.

The environment will be ready for use. You can now run PyPSA-Earth in the container.
