FROM condaforge/mambaforge

RUN conda update -n base conda
RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba

WORKDIR /pypsa-earth

COPY ./envs ./envs

RUN conda env create -n pypsa-earth --file envs/environment.yaml

RUN rm -r envs

RUN echo "conda activate pypsa-earth" > ~/.bashrc
ENV PATH /opt/conda/envs/pypsa-earth/bin:$PATH

ENTRYPOINT ["tail", "-f", "/dev/null"]
