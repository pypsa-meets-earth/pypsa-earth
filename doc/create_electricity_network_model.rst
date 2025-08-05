.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

##########################################
Building electricity networks
##########################################

The simplification ``snakemake`` rules prepare **approximations** of the full model, for which it is computationally viable to co-optimize generation, storage and transmission capacities.

- :mod:`base_network` builds and stores the base network with all buses, HVAC lines and HVDC links
- :mod:`add_electricity` adds the generators and demand to the network model
- :mod:`simplify_network` transforms the transmission grid to a 380 kV only equivalent network
- :mod:`cluster_network` uses a clustering approach (e.g. `k-means <https://en.wikipedia.org/wiki/K-means_clustering>`_ )
  to partition the network into a given number of zones and then reduce the network to a representation with one bus per zone.
- :mod:`add_extra_components` add extra components to the model, such as storage
- :mod:`prepare_network` introduces optional constraints and requirements in the modelling,
  such as CO2 emissions, security margins, etc.

The simplification and clustering steps are described in detail in the paper

- Jonas Hörsch and Tom Brown. `The role of spatial scale in joint optimisations of generation and transmission for European highly renewable scenarios <https://arxiv.org/abs/1705.07617>`_), *14th International Conference on the European Energy Market*, 2017. `arXiv:1705.07617 <https://arxiv.org/abs/1705.07617>`_, `doi:10.1109/EEM.2017.7982024 <https://doi.org/10.1109/EEM.2017.7982024>`_.

Index:

.. toctree::
   :caption: Overview

   create_network_model/base_network
   create_network_model/add_electricity
   create_network_model/simplify_network
   create_network_model/cluster_network
   create_network_model/add_extra_components
   create_network_model/prepare_network
