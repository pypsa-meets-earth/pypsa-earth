# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: CC-BY-4.0
FROM condaforge/mambaforge

RUN conda update -n base conda
RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba

RUN apt-get update && apt-get install -y bash git && apt-get install gcc -y

WORKDIR /pypsa-earth

COPY ./envs ./temp

RUN conda env create -n pypsa-earth -f temp/linux-64.lock.yaml

RUN conda init bash

RUN touch ~/.bashrc && echo "conda activate pypsa-earth" >> ~/.bashrc

SHELL ["/bin/bash", "--login", "-c"]

ENV PATH /opt/conda/envs/pypsa-earth/bin:$PATH

RUN conda install conda-forge::openjdk -y

RUN rm -r temp

RUN conda clean -afy && \
    rm -rf /tmp/*

CMD ["bash"]
