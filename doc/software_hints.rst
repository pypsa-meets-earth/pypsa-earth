.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _software_hints:


##########################################
Software hints
##########################################

Tools we are using
=====================

The complete list of software needed for work with PyPSA Earth is listed below.

Python
-----------------------------------

`Python 3 <https://www.python.org/>`_ is used as our main programming language, thus its knowledge is mandatory.

.. TODO Add installation instructions

To refresh your knowledge, you may use plenty of online courses free-of-charge, e.g. `CSDojo playlist <https://www.youtube.com/c/CSDojo/playlists>`_. Useful content to watch refer to numpy and pandas.

Conda
-----------------------------------

In order to use packages in python, it is highly recommended to use a `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html>`_ package manager, such as `Anaconda <https://docs.anaconda.com/>`_ 

You may check if ``conda`` is already installed on your system with `conda -V`. In case it's not we recommend to install `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_, which is a mini version of `Anaconda <https://www.anaconda.com/>`_ and includes only ``conda`` and its dependencies. 

For instructions for your operating system follow the ``conda`` `installation guide <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_.

We recommend to install `mamba <https://github.com/QuantStack/mamba>`_ in `conda` base environment as that makes dependencies management much faster. 

Useful hints on installing Anaconda and setting up Ubuntu terminal with VSCode in Windows is available `here <https://gist.github.com/kauffmanes/5e74916617f9993bc3479f401dfec7da>`_.

There are many things which could go wrong with conda. `This article <https://towardsdatascience.com/conda-essential-concepts-and-tricks-e478ed53b5b>`_ provides you a crystal clear explanation of conda (**excellent read**).
 
Git
-----------------------------------

`Git <https://git-scm.com/>`__ is a free open source system aimed at tracking changes in the code development and enable to coordinate the parallel software development between many developers.

It is highly recommended to `learn the git basics <https://git-scm.com/doc>`_.

.. TODO Add Git tutorials


.. Not sure if it's needed 
.. Java
.. ----------------------

.. `Java <https://www.oracle.com/java/technologies/downloads/>` is needed for using `powerplantmatching` package. To have a better user experience, please install the redistribution from the website according to your operating system.

Solvers
-----------------------------------

Useful instructions of installation and setup:

* `HiGHS solver <https://github.com/PyPSA/PyPSA/blob/633669d3f940ea256fb0a2313c7a499cbe0122a5/pypsa/linopt.py#L608-L632>`_;

* `Gurobi solver <https://www.youtube.com/watch?v=yNmeG6Wom1o>`_.
 
Integrated Development Environment
-----------------------------------

In order to write and debug python code, you need an Integrated Development Environment (IDE) that is a software used to write code. We recommend `Visual Studio Code <https://code.visualstudio.com/>`_, which is freely available online and provides an easy to use interface with Git. Obviously, any alternatives like `PyCharm <https://www.jetbrains.com/pycharm/>`_ or `Sublime <https://www.sublimetext.com/>`_ will work as well.


System requirements
===================

Building the model with the scripts in this repository runs on a normal computer with 8-16 GB RAM. Depending of the region of interest, different amounts of storage (HHD/SSD) are required. Model of Africa requires about 40 Gb, the world 250 Gb, a single country between 1-10 Gb. We have also prepared a tutorial which should be below 10Gb.
