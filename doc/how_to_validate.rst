..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _tutorial::

##########################################
How to validate?
##########################################

.. TODO add a list of actions needed to do the validation

To validate the data obtained with PyPSA-Earth, we recommend to go through the procedure here detailed. An exampled of the validation procedure is available in the `Nigeria validation <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/validation_nigeria.ipynb>`_ notebook. Public information on the power system of Nigeria are compared to those obtained from the PyPSA-Earth model.

Simulation procedure
===================================

It may be recommended to check the following quantities the validation:

#. inputs used by the model:

    #. network characteristics;

    #. substations;

    #. installed generation by type;

#. outputs of the simulation:

    #. demand;

    #. energy mix.

Data sources
===================================
 
Data availability for many parts of the world is still quite limited. Usually the best sources to compare with are regional data hubs. There is also a collection of harmonized datasets curated by the international organisations. A non-exhaustive list of helpful sources:

* `World Bank <https://energydata.info/>`_;

* International Renewable Energy Agency `IRENA <https://pxweb.irena.org/pxweb/en/IRENASTAT/IRENASTAT__Power%20Capacity%20and%20Generation/ELECCAP_2022_cycle2.px/>`_;

* International Energy Agency `IEA <https://www.iea.org/data-and-statistics>`_;

* `BP <https://www.bp.com/en/global/corporate/energy-economics/statistical-review-of-world-energy.html>`_ Statistical Review of World Energy;

* International Energy Agency `IEA <https://www.iea.org/data-and-statistics>`_;

* `Ember <https://ember-climate.org/data/data-explorer/>`_ Data Explorer.


Advanced validation examples
===================================

The following validation notebooks are worth a look when validating your energy model:

1. A detailed `network validation <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/network_validation.ipynb>`_.
 
2. Analys of `the installed capacity <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/capacity_validation.ipynb>`_ for the considered area. 

3. Validation of `the power demand <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/demand_validation.ipynb>`_ values and profile.

4. Validation of `hydro <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/hydro_generation_validation.ipynb>`_, `solar and wind <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/renewable_potential_validation.ipynb>`_ potentials.
