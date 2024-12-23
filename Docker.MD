<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: AGPL-3.0-or-later
-->
# PyPSA-Earth Development Environment Setup with Docker

This guide provides an alternative way to set up a development environment for **PyPSA-Earth** using Docker containers. This method is particularly useful for users who are not familiar with programming or Python, or who prefer not to install Python directly on their local machine. Using Docker simplifies setup by creating a self-contained environment for PyPSA-Earth.

## Prerequisites

### 1. Install Docker

Ensure Docker is installed on your system. Follow the instructions for your operating system:

- **Windows**: [Docker for Windows](https://docs.docker.com/desktop/install/windows-install/)
- **Linux**: [Docker for Linux](https://docs.docker.com/desktop/install/linux/)
- **MacOS**: [Docker for MacOS](https://docs.docker.com/desktop/install/mac-install/)

### 2. Install Visual Studio Code(VSC)

You can use the link [here](https://code.visualstudio.com/download) to install Visual Studio Code on your operating system. Ensure to select the most compatible file for your operating system.

### 3. Install GitHub Desktop

You will also need GitHub Desktop to clone the PyPSA-Earth repository. Install GitHub Desktop for your operating system from [here](https://desktop.github.com/download/).

## Steps to Set Up PyPSA-Earth with Docker

### Step 1: Clone the Repository

1. Open **GitHub Desktop**.
2. Go to **File** > **Clone Repository**.
3. Paste the following URL in the URL field:

    ```bash
    https://github.com/pypsa-meets-earth/pypsa-earth.git
    ```

4. Click on **Clone**.
5. Choose the location where you want to save the repository on your local machine.
6. After cloning, click on **Current Branch: main** and select `devContainers`.
7. Click on **Open in Visual Studio Code**.

### Step 2: Rebuild and Open in Container

1. Open the cloned repository in **VSCode**.
2. Click on the icon located at the bottom left corner of the VSCode window.
3. Select **Reopen in Container**.
4. Wait for the container to build and for the repository to open in the container.

Once these steps are completed, your development environment will be ready, and you can start using **PyPSA-Earth** in the Docker container.

---

You are now all set up! You can use the development environment to explore PyPSA-Earth and make modifications as needed within the Docker container.

You can start running the tutorial [here](https://pypsa-earth.readthedocs.io/en/latest/short_tutorial.html)
