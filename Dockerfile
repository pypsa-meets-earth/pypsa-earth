FROM condaforge/mambaforge

RUN conda update -n base conda
RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba

RUN apt-get update && apt-get install -y bash git

WORKDIR /pypsa-earth

COPY ./envs ./temp

RUN conda env create -n pypsa-earth -f temp/environment.yaml

RUN conda init bash

RUN touch ~/.bashrc && echo "conda activate pypsa-earth" >> ~/.bashrc

SHELL ["/bin/bash", "--login", "-c"]

ENV PATH /opt/conda/envs/pypsa-earth/bin:$PATH

RUN conda install conda-forge::openjdk -y

RUN rm -r temp

RUN conda clean -afy && \
    rm -rf /tmp/*

CMD ["bash"]
