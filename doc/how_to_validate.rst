..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _tutorial::

##########################################
How to validate?
##########################################

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