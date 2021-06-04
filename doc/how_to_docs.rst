..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Africa authors

  SPDX-License-Identifier: CC-BY-4.0

.. _how_to_docs:

##########################################
How to docs?
##########################################

We add the code documentation along the way. You might think that cost a lot of time and is not efficient - but thats not really true anymore! Documenting with great tools makes life much easier for YOU and YOUR COLLABORATORS and speed up the overall process. Using `Readthedocs <https://docs.readthedocs.io/en/stable/intro/getting-started-with-sphinx.html>`_ and its add on `sphinx.ext.autodoc  <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_ we document in our code scripts which then will automatically generate the documentation you might see here. 

Thank you Eric Holscher & team for your wonderful *Readthedocs* open source project. You can find an emotional speak by Eric `here <https://www.youtube.com/watch?v=U6ueKExLzSY>`_. 

Structure and Syntax example 
=============================

The documentation is fully stored in our `doc` folder. You might note that most files are used as *reStructuredText* file or .rst which is very popular by coders when documentation matters. You could differantiate between to elements:

1. Non-automated doc elements. They simply make the text appealing. Example in `installation.rst` in our doc folder. To write these requires some knowledge on writing the text which is quite easy to learn having this `**cheat sheet** <https://github.com/DevDungeon/reStructuredText-Documentation-Reference#syntax-examples>`_ in close reach.
2. Automated doc elements using automodule\class or similar. What they do is basically to link the code script with the doc texts, for instance, **compare the python script** `osm_pbf_power_data_extractor.py <https://github.com/pz-max/pypsa_meets_africa/blob/main/data_exploration/osm_pbf_power_data_extractor.py>`_ Example `api_reference.rst` with the documentation. To write these kind of automation, get inspiration from our or the `pypsa doc for the api_reference <https://pypsa.readthedocs.io/en/latest/api_reference.html>`_ documentation. Further, to help understanding how things work the official `Sphinx documentation <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_ might help too. 

We found three important files/file groups for the documentation:
1. `index.rst`. This is the starter page of the *Readthedocs* page and generates also the sidebar with the outline (see `pypsa <https://pypsa.readthedocs.io/en/latest/index.html>`_).
2. `The .py script` with the acutal code documentation.
3. `conf.py` with some useful information if something does not work.

Please, if you have problems with the documentation create an issue and let us know

