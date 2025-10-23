.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _customization_validation:

###################
4. Model Validation
###################

To validate the data obtained with PyPSA-Earth, we recommend to go through the procedure here detailed. An exampled of the validation procedure is available in the `Nigeria validation <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/validation_nigeria.ipynb>`_ notebook. Public information on the power system of Nigeria are compared to those obtained from the PyPSA-Earth model.

Data workflow
^^^^^^^^^^^^^^

PyPSA-Earth is designed to provide an optimisation model workflow with all necessary data for every country of the world. An architecture of the data workflow looks like follows.

#. Inputs used by the workflow:

    #. Aggregated power and energy demand by sectors

    #. Installed generation capacity including it's spatial distribution

    #. Energy-relevant weather parameters, land usage and orthography parameters

#. Workflow outputs fed into the optimisation model:

    #. Network topology and transmission capacity

    #. Spatially distributed power and energy demands

    #. Available renewable potential

#. Outputs of the optimisation run:

    #. Power and energy mix

    #. Dispatch time-series

    #. Capital and operational costs

It may be recommended to check the major quantities in a validation procedure both for the workflow and the optimisation.

Where to look for reference data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Data availability for many parts of the world is still quite limited. Usually the best sources to compare with are regional data hubs. There is also a collection of harmonized datasets curated by the international organisations. A non-exhaustive list of helpful sources:

* `Ember <https://ember-climate.org/data/data-explorer/>`_ Data Explorer

* `World Bank <https://energydata.info/>`_

* International Renewable Energy Agency `IRENA <https://pxweb.irena.org/pxweb/en/IRENASTAT/>`_

* International Energy Agency `IEA <https://www.iea.org/data-and-statistics>`_

* A global collection of power grid charts `Awesome Electrical Grid Mapping <https://github.com/open-energy-transition/Awesome-Electrical-Grid-Mapping>`_

* `BP <https://www.bp.com/en/global/corporate/energy-economics/statistical-review-of-world-energy.html>`_ Statistical Review of World Energy

How to tailor the inputs
^^^^^^^^^^^^^^^^^^^^^^^^
1. Adjust configuration parameters
2. Use custom data
3. Improve data quality in the source

Contribution guide

`Visualization toolkit <https://github.com/ben10dynartio/osm-power-grid-map-analysis>`__ to quickly get insights into OSM-extracted power grid datausing maps and graph analysis

Powerplantmatching

GEM

atlite

ERA5

Advanced validation examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following validation notebooks are worth a look when validating your energy model:

1. A detailed `network validation <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/network_validation.ipynb>`_.

2. Analysis of `the installed capacity <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/capacity_validation.ipynb>`_ for the considered area.

3. Validation of `the power demand <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/demand_validation.ipynb>`_ values and profile.

4. Validation of `hydro <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/hydro_generation_validation.ipynb>`_, `solar and wind <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/renewable_potential_validation.ipynb>`_ potentials.
