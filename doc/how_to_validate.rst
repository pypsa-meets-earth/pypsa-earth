..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _tutorial::

##########################################
How to validate?
##########################################

Data workflow is PyPSA-Earth allows to extract all the data needed to build your energy model from global open data sources. However, data availability for many parts of the world is still quite limited and it's essential to validate the data included into the model. Usually the best sources to compare with are regional data hubs along with harmonized datasets of international organisations such as 
- `World Bank <https://energydata.info/>`_;
- International Renewable Energy Agency `IRENA <https://pxweb.irena.org/pxweb/en/IRENASTAT/IRENASTAT__Power%20Capacity%20and%20Generation/ELECCAP_2022_cycle2.px/>`_;
- International Energy Agency `IEA <https://www.iea.org/data-and-statistics>`_
- `BP Statistical Review of World Energy <https://www.bp.com/en/global/corporate/energy-economics/statistical-review-of-world-energy.html>`_;
- International Energy Agency `IEA <https://www.iea.org/data-and-statistics>`_

The following validation points are worth keeping in mind when validating your energy model:

1. Check the `power grid <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/network_validation.ipynb>`_:
    - overall lines length;
    - general grid topology;
    - ensure that the general structure of the grid model is appropriate, playing with `tol` values and augmentation options if needed.
 
2. Compare the `installed capacity <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/capacity_validation.ipynb>`_ values 

3. Validate the `power demand <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/demand_validation.ipynb>`_ values and profile.

4. Check that `hydro <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/hydro_generation_validation.ipynb>`_, `solar and wind <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/renewable_potential_validation.ipynb>`_ potentials have reasonable values

5. Simulate the actual `energy mix <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/validation_nigeria.ipynb>`_. Look for detailed explanations in https://arxiv.org/abs/2209.04663, section 5.1.

Data availability and quality usually is the biggest concern. Some useful hints on the real-world validation example can be found in the `Nigeria validation <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/validation_nigeria.ipynb>`_ notebook.