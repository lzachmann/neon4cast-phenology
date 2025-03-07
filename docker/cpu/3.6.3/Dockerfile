FROM rocker/tidyverse:3.6.3

ENV WORKON_HOME /opt/virtualenvs
ENV PYTHON_VENV_PATH $WORKON_HOME/r-tensorflow

## Set up a user modifyable python3 environment
RUN apt-get update && apt-get install -y --no-install-recommends \
        libpython3-dev \
        python3-venv && \
    rm -rf /var/lib/apt/lists/*

RUN python3 -m venv ${PYTHON_VENV_PATH}

RUN chown -R rstudio:rstudio ${WORKON_HOME}
ENV PATH ${PYTHON_VENV_PATH}/bin:${PATH}
## And set ENV for R! It doesn't read from the environment...
RUN echo "PATH=${PATH}" >> /usr/local/lib/R/etc/Renviron && \
    echo "WORKON_HOME=${WORKON_HOME}" >> /usr/local/lib/R/etc/Renviron && \
    echo "RETICULATE_PYTHON_ENV=${PYTHON_VENV_PATH}" >> /usr/local/lib/R/etc/Renviron

## Because reticulate hardwires these PATHs...
RUN ln -s ${PYTHON_VENV_PATH}/bin/pip /usr/local/bin/pip && \
    ln -s ${PYTHON_VENV_PATH}/bin/virtualenv /usr/local/bin/virtualenv

## install as user to avoid venv issues later
USER rstudio
RUN pip3 install \
    h5py==2.9.0 \
    pyyaml==3.13 \
    requests==2.21.0 \
    Pillow==5.4.1 \
    tensorflow==1.14.0 \
    tensorflow-probability==0.7.0 \
    keras==2.2.4 \
    --no-cache-dir
USER root
RUN install2.r reticulate tensorflow keras

# Install ncdf4 package.
RUN apt-get update && apt-get install -y libnetcdf-dev && \
    wget http://cirrus.ucsd.edu/~pierce/ncdf/ncdf4_1.13.tar.gz && \
    R CMD INSTALL ncdf4_1.13.tar.gz --configure-args="--with-nc-config=/usr/bin/nc-config" && \
    rm ncdf4_1.13.tar.gz

## greta and other greenwave dependencies
RUN install2.r --error \
    aws.s3 \
    abind \
    moments \
    mvtnorm \
    bayesplot \
    DiagrammeR \
    ggthemes \
    HDInterval \
    NCmisc \
    txtplot \
    webshot && \
  Rscript -e "update.packages(ask = FALSE)" && \
  R -e "devtools::install_github('greta-dev/greta')"
