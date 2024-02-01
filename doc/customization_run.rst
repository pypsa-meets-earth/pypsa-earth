.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _customization_run:

####################
2. Running the model
####################

To solve full optimization problem, it is important to pick a solver. For instance, this tutorial uses the open-source solver glpk and does not rely
on the commercial solvers such as Gurobi or CPLEX (for which free academic licenses are available).

.. code:: yaml

    solver:
        name: glpk

.. note::

    ``glpk`` can solve the network with low temporal and spacial resolution. To make a full model run, it is adviced to use ``CPLEX``, ``Gurobi``, or open-source `HIGHs <https://highs.dev/>`__.



