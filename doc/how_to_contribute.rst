..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _how_to_contribute:

##########################################
How to contribute
##########################################

Contributions are welcome, and they are greatly appreciated! 
Every little bit helps, and you always earn credits.

You can contribute on the code side in many ways:

- submit feedback,
- add new features,
- report bugs, 
- fix bugs, 
- implement a new cluster/cloud computation backend,
- write documentation


Code
====

For every code contribution you should run `pre-commit <https://pre-commit.com/index.html>`_.
This will lint, format and check your code contributions against our guidelines
(e.g. we use `Black <https://github.com/psf/black>`_ as code style
and aim for `REUSE compliance <https://reuse.software/>`_):

1. Installation ``conda install -c conda-forge pre-commit`` or ``pip install pre-commit``
2. Usage:
  * To automatically activate ``pre-commit`` on every ``git commit``: Run ``pre-commit install``
  * To manually run it: ``pre-commit run --all``

Add a new test if you want to contribute new functionality to the config.
We perform currently *multiple* integration tests which means various workflows need to work.
All test configs are build by updating the `config.tutorial.yaml` with the configs in `pypysa-earth/test/*.yaml`.

  * You can test your contribution locally with `snakemake --cores 4 run_tests`. This will build test configs and executes them.
  * Run `snakemake -j1 build_test_configs` to build and analyse locally the test configs.

To contribute a test:

1. Provide a new test in `test/<new test>.yaml`, or adjust one of the existing ones. **config.standard** updates the config.tutorial.yaml, all other existing or new tests, e.g. custom or landlock, are built from the *standard test*. 
2. Add a new test config path to the *rule build_all_test* in the *Snakefile*.
3. If your functionality should be tested in the CI for every pull request, add a respective code in `.github/workflows/ci*.yaml`


No-Code
========
Instead of contributing code there are alternatives to support the PyPSA-Earth goals.
You can fund projects, supervise people, support us with outreach activities or events.
Check out our `website <https://pypsa-meets-earth.github.io>`_ for more details.


Join us and get involved
========================

Any person/ group is welcome to join us. Be it research leader, researcher, undergraduate, or industry professional.
A simple way to explore opportunities for collaboration is to join our meetings. All of them are **OPEN**.

- `List of meetings and times <https://github.com/pypsa-meets-earth/pypsa-earth#get-involved>`_

- **Discord**
  
  - Chat with the community, team up on features, exchange with developers, code in voice channels
  - `Discord invitation link <https://discord.gg/AnuJBk23FU>`_
