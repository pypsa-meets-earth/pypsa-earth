.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _data:

Rule ``retrieve_databundle_light``
==================================

Not all data dependencies are shipped with the git repository, since git is not suited for handling large changing files. Instead we provide separate data bundles which can be obtained using the ``retrieve_databundle_light`` rule when `retrieve_databundle` flag in the configuration file is on. If that is the case, `retrieve_databundle_light` rule is included into the workflow. The common data needed to run the model will be loaded corresponding to settings of the `config_default.yaml` or `config_tutorial.yaml` depending on the `tutorial` flag. 

.. automodule:: retrieve_databundle
    :noindex:
