.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _learning_materials:

Learning materials
===================================

PyPSA meets Earth builds on top of several open-source packages, which are here recalled together with recommended sources to learn them from scratch.

.. _data_science_basics:

YouTube Tutorials
---------------------

We host **training** and **presentation** videos about PyPSA ecosystem tools made for users, developers and interested people.

- `Training playlist <https://www.youtube.com/playlist?list=PLz1RcE1NdLlEzCT8G_v-sJKOeHOhZvw5R>`_
- `Presentation playlist <https://www.youtube.com/playlist?list=PLz1RcE1NdLlHpsyopZUKDkeyAuQotfOE7>`_
- `Youtube channel <https://www.youtube.com/@pypsameetsearth7475/featured>`_

External community resources

- `Energy System Modelling Playlist <https://www.youtube.com/playlist?list=PLaXo6_ws0ZKJUYyUu-e6HRUcjuZnen7pA>`_ by Anton Achhammer: An in-depth series exploring advanced energy system modeling techniques and insights from an experienced researcher in the field.


PyPSA Introduction (essential)
-------------------------------

- `Watch <https://www.youtube.com/watch?v=ty47YU1_eeQ>`_ how PyPSA-Eur is designed (1h)
- `Watch <https://www.youtube.com/watch?v=mAwhQnNRIvs>`_ and put your hands on to make PyPSA-Eur work on your computer (1-3h)
- While watching these PyPSA videos always have a look into the excellent `PyPSA-Eur documentation <https://pypsa-eur.readthedocs.io/en/latest/index.html>`_
- To see what data we can extract we work usually closely with the `basic PyPSA documentation <https://pypsa.readthedocs.io/en/latest/components.html>`_
- Follow the `PyPSA-Earth hackathon <https://github.com/pypsa-meets-earth/documentation#1-hackathon-material>`_. Check the videos from this `Youtube series <https://www.youtube.com/watch?v=2nakM174CME&list=PLz1RcE1NdLlEzCT8G_v-sJKOeHOhZvw5R&index=7>`_ since it explains and applies partially the hackathon material.

Data science basics (essential)
--------------------------------

- Fabian Neumann just shared with the world the possibly best training material for `"Data Science for Energy System Modelling" <https://fneum.github.io/data-science-for-esm/intro.html>`_. This is a free multi-week course preparing you for all you need for PyPSA-Earth.
- Refresh your Python knowledge by watching `CSDojo's playlist <https://www.youtube.com/c/CSDojo/playlists>`_. His content is excellent as introduction. You will learn in effective short videos the python basics such as variables If/else statements, functions, lists, for loops, while loops, dictionaries, classes and objects, boolean, list comprehensions, sets - put your hands on and write some test scripts as the video suggests. (~3h)
- Familiarize yourself with numpy and panda dataframes.  In the Python-based PyPSA tool, we do not work with Excel. Powerful panda dataframes are our friends. `Here <https://www.coursera.org/learn/python-data-analysis>`__ is an extensive 30h course that provides a great introduction if this is completely unfamiliar to you.
- `Introduction to Unix-shell <https://swcarpentry.github.io/shell-novice/>`_ - "Use of the shell is fundamental to a wide range of advanced computing tasks, including high-performance computing and automated workflow. These lessons will introduce you to this powerful tool." (optional 4h, to become a pro)


Git and GitHub
--------------

For code collaboration we use GitHub. Which is a common source control tool that is a very popular collaborative code development tool. Here some notes if you are not already familiar with it:

- Git and GitHub is not the same. Usually, you work with git on your computer (offline) to push changes to GitHub (online).
- `Here <https://www.youtube.com/watch?v=8JJ101D3knE>`__ a great intro which we recommend
- Learning by doing. Maybe one of the best ways to learn is to puts your hands on open a GitHub repository and upload/change/reverse files from your local computer on some dummy scripts.
- This `cheatsheet <https://www.atlassian.com/git/tutorials/atlassian-git-cheatsheet>`_ might help using the Git commands


Snakemake and advanced changes
------------------------------

Snakemake is our brain in PyPSA.
It automates many tasks & keeps the code structure clean.
Therefore, it is quite useful to learn if your task is to integrate features into PyPSA.
We can recommend:

- `snakemake basic and advanced tutorial here <https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html>`__ (takes max 3-5h and makes a lot of fun).
- Explore how PyPSA uses snakemake in the Snakefile and scripts - the GitHub search function is your best friend to find quickly what interests you.


Code environment
----------------

We can recommend setting up VScode from Microsoft. Add some extension if you like as described in `this video <https://www.youtube.com/watch?v=0fROnrISdZU>`_. For instance GitHub, Gitlense, and maybe some others.

*Note*: if you decide to use Visual Studio Code, check out the tutorial about how to use `Git <https://code.visualstudio.com/docs/editor/versioncontrol#_git-support>`_ and `Github <https://code.visualstudio.com/docs/editor/github>`_  in Visual Studio Code
