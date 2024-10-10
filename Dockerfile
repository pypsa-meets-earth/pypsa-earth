FROM condaforge/mambaforge

RUN conda update -n base conda
RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba

WORKDIR /pypsa-earth

COPY ./envs ./envs

RUN conda env create --file envs/environment.yaml

RUN rm -r envs

RUN echo "source activate pypsa-earth" > ~/.bashrc
ENV PATH /opt/conda/envs/pypsa-earth/bin:$PATH

ENTRYPOINT ["tail", "-f", "/dev/null"]
