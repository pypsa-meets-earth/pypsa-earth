.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _how_to_contribute:

##########################################
How to Сontribute
##########################################

Contributions are welcome, and they are greatly appreciated!
Every little bit helps, and you always earn credits.

You can contribute on the code side in many ways:

* submit feedback,
* add new features,
* report bugs,
* fix bugs,
* implement a new cluster/cloud computation backend,
* write documentation


Code
====

Linting and pre-commit
----------------------
For every code contribution you should run `pre-commit <https://pre-commit.com/index.html>`_.
This will lint, format and check your code contributions against our guidelines
(e.g. we use `Black <https://github.com/psf/black>`_ as code style
and aim for `REUSE compliance <https://reuse.software/>`_):

1. Installation ``conda install -c conda-forge pre-commit`` or ``pip install pre-commit``
2. Usage:
  * To automatically activate ``pre-commit`` on every ``git commit``: Run ``pre-commit install``
  * To manually run it: ``pre-commit run --all``

Testing
-------
Add a new test if you want to contribute new functionality to the config.
We perform currently *multiple* integration tests which means various workflows need to work.
All test configs are build by updating the ``config.tutorial.yaml`` with the configs in ``pypysa-earth/test/*.yaml``.

  * You can test your contribution locally with ``make test``.
  * See the Makefile for further information which configurations are tested.

To contribute a test:

1. Provide a new test in ``test/<new test>.yaml``, or adjust one of the existing ones. These tests update the config.tutorial.yaml to test other options e.g. landlock countries.
2. Add a new test config path to the ``rule build_all_test`` in the ``Snakefile``.
3. If your functionality should be tested in the CI for every pull request, add a respective code in ``.github/workflows/ci-linux.yaml``. We test all functionalities only for Linux while providing a general test for windows and mac.

Performance-profiling
---------------------
Performance profiling is important to understand bottlenecks and
the accordingly optimize the speed in PyPSA-Earth. We use the Python built-in
`cProfiler`, custom decorators on single functions and analysis tools
like `snakeviz <https://jiffyclub.github.io/snakeviz/>`_. See a detailed example
in `this discussion #557 <https://github.com/pypsa-meets-earth/pypsa-earth/discussions/557>`_.


Documentation
============
How to docs?
-------------

We add the code documentation along the way.
You might think that cost a lot of time and is not efficient - but that's not really true anymore!
Documenting with great tools makes life much easier for YOU and YOUR COLLABORATORS and speed up the overall process.
Using `Readthedocs <https://docs.readthedocs.io/en/stable/intro/getting-started-with-sphinx.html>`_ and its add-on
`sphinx.ext.autodoc  <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_ we document in our
code scripts which then will automatically generate the documentation you might see here.

Thank you Eric Holscher & team for your wonderful *Readthedocs* open source project.
You can find an emotional speech by Eric `here <https://www.youtube.com/watch?v=U6ueKExLzSY>`_.

Structure and Syntax example
-----------------------------

The documentation is fully stored in our `doc` folder. You might note that most files are used as *reStructuredText* file or .rst which is very popular by coders when documentation matters.

You could differentiate between to elements:

1. Non-automated doc elements. They simply make the text appealing. Example in `installation.rst` in our doc folder. To write these requires some knowledge on writing the text which is quite easy to learn having this `**cheat sheet** <https://github.com/DevDungeon/reStructuredText-Documentation-Reference#syntax-examples>`_ in close reach.
2. Automated doc elements using automodule\class or similar. What they do is basically to link the code script with the doc texts, for instance, **compare the python script** `add_electricity.py <https://github.com/pypsa-meets-earth/pypsa-earth/blob/main/scripts/add_electricity.py>`_ Example `api_reference.rst` with the documentation. To write these kind of automation, get inspiration from our or the `pypsa doc for the api_reference <https://pypsa.readthedocs.io/en/latest/api_reference.html>`_ documentation. Further, to help understanding how things work the official `Sphinx documentation <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_ might help too.

We found three important files/file groups for the documentation:

1. `index.rst`. This is the starter page of the *Readthedocs* page and generates also the sidebar with the outline (see `pypsa <https://pypsa.readthedocs.io/en/latest/index.html>`_).
2. `The .py script` with the actual code documentation.
3. `conf.py` with some useful information if something does not work.

The images for documentation should be placed into `documentation <https://github.com/pypsa-meets-earth/documentation>`_ repository to the folder "doc/img". The content of the folder "documentation/doc/img/" is copied into "pypsa-earth/doc/img/" during building PyPSA-Earth documentation.

Please, if you have problems with the documentation create an issue and let us know


To create the documentation locally, you need `Sphinx <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_ . It can be installed using specifications
from `doc/requirements.txt`. First, we recommend creating a fresh conda environment and activate it:

.. code:: bash

    .../pypsa-earth % conda create --name pypsa-earth-docs python

    .../pypsa-earth % conda activate pypsa-earth-docs

Next, install the packages specified in `doc/requirements.txt` using `pip`:

.. code:: bash

    .../pypsa-earth % pip install -r doc/requirements.txt


Once installation is completed, the following commands allow you to create the documentation locally:

.. code:: bash

    .../pypsa-earth (pypsa-earth-docs) % cd doc

    .../pypsa-earth/doc (pypsa-earth-docs) % make html

This will create html files in `pypsa-earth/doc/_build/html`.
VScode provides a so called Liveserver extension such that the html file can be opened locally on your computer.

.. warning::

    Windows users might face some challenges when building the documentation locally using `make`. A workaround can be found, but might be time consuming. For instance:

    1. If using Windows PowerShell, one might need to replace the command `make html` above by `./make html`. For more details on what is going on, see `this post <https://stackoverflow.com/questions/65471557/make-html-not-working-for-sphinx-documentation-in-windows-10>`_ on Stack Overflow.

.. note::

    The documentation is built automatically by the CI for every pull request. The documentation is hosted on `ReadTheDocs <https://pypsa-earth.readthedocs.io/en/latest/>`_.
    For more information on our documentation infrastructure and syntax tips, see `this page <https://pypsa-earth.readthedocs.io/en/latest/how_to_docs.html>`_.


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
