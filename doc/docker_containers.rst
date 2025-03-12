.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _docker_containers:

Installation with Docker
===============================================

This is an alternative way to create a development environment for PyPSA-Earth. This method is useful for users who are not familiar with programming or Python, or who do not want to install Python on their local machine. It uses Docker containers to create a development environment for PyPSA-Earth.

This section provides a step-by-step guide on how to set up and use Docker containers to run PyPSA-Earth.

Steps:

1. Install Docker: Follow the instructions for your operating system:

* `Windows <https://docs.docker.com/desktop/install/windows-install/>`_
* `Linux <https://docs.docker.com/desktop/install/linux/>`_
* `MacOS <https://docs.docker.com/desktop/install/mac-install/>`_

    Ensure Docker is installed on your system.

2. You can use the link `here <https://code.visualstudio.com/download>`_ to install Visual Studio Code on your operating system. Ensure to select the most compatible file for your operating system.

3. Install GitHub Desktop for your OS `here <https://desktop.github.com/download/>`_.

4. Clone the repository:
    * Open GitHub Desktop.
    * Click on "File" in the top left corner.
    * Click on "Clone Repository".
    * Paste the following URL in the URL field:

    .. code:: bash

        https://github.com/pypsa-meets-earth/pypsa-earth.git

    * Click on "Clone".
    * Choose the location where you want to save the repository.
    * Click on "Current Branch: main" and select `devContainers`.
    * Click on "Open in Visual Studio Code".

    The repository will be cloned to your local machine.

5. Rebuild and open in a container:
    * Open the repository in VSCode.
    * Click on the icon in the far bottom left corner of the VSCode window.
    * Click on "Reopen in Container".
    * Wait for the container to build and open the repository in the container.

    The environment will be ready for use. You can now run PyPSA-Earth in the container.
