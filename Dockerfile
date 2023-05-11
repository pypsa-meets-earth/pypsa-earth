FROM --platform=linux/x86_64 mambaorg/micromamba

WORKDIR /pypsa-earth

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/environment.yaml /tmp/env.yaml

USER root

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git

COPY gurobi.lic /opt/gurobi/gurobi.lic
ENV GRB_LICENSE_FILE=/opt/gurobi/gurobi.lic

RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python will not be found)



RUN python -c 'import uuid; print(uuid.uuid4())' > /tmp/my_uuid

RUN python -c "import pypsa"

COPY . /pypsa-earth
