FROM condaforge/mambaforge

ENV INSTALL_DIR="/opt"

RUN apt-get update; export DEBIAN_FRONTEND=noninteractive; apt-get -y --no-install-recommends install \
   wget \
   unzip \
   build-essential \
   gcc-multilib \
   multiqc \
   netcat-openbsd && \
   export DEBIAN_FRONTEND=noninteractive; apt-get -y install openssh-client openssh-sftp-server openssh-server

# Install miniconda
# ENV CONDA_DIR /opt/conda
# RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
#    /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH
# Install the conda environments
COPY * /
RUN ls -al
RUN sh installEnvs.sh

#RUN source ~/opt/conda/etc/profile.d/conda.sh
USER 1000
