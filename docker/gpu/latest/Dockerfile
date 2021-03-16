FROM tensorflow/tensorflow:1.14.0-gpu-py3-jupyter

LABEL org.label-schema.license="CSP-1.0" \
    org.label-schema.vcs-url="https://gitlab.com/tonychang" \
    org.label-schema.vendor="Conservation Science Partners" \
    maintainer="Tony Chang <tony@csp-inc.org>"

#install R
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        apt-utils \
        apt-transport-https \
        ca-certificates \
        cron \
	curl \
        ed \
	gdebi-core \
	git \
        gfortran \
        gsfonts \
        gnupg2 \
        less \
        libapparmor1 \
        libblas-dev \
        libbz2-dev \
        libclang-dev \
        libcurl4-openssl-dev \
        libegl1-mesa \
        libfreetype6-dev \
        libhdf5-serial-dev \
        libicu-dev \
        libpcre2-dev \
        libpcre3-dev \
        liblapack-dev \
        liblzma-dev \
        libssl1.0.0 \
        libssl-dev \
	libuser \
	libuser1-dev \
        libxml2-dev\
        libxcursor-dev \
        libxkbcommon-x11-0 \
        libzmq3-dev \
        locales \
	pkg-config \
        psmisc \
	rrdtool \
        software-properties-common \
	sudo \
        texinfo \
        texlive-extra-utils \
        texlive-fonts-recommended \
        texlive-fonts-extra \
        texlive-latex-recommended \
        unzip \
        vim-tiny \
        wget \
        && add-apt-repository --enable-source --yes "ppa:marutter/rrutter3.5" \
	&& add-apt-repository --enable-source --yes "ppa:marutter/c2d4u3.5"

#install nodejs
RUN curl -sL https://deb.nodesource.com/setup_10.x | bash
RUN apt-get install -y nodejs default-jdk default-jre

# Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8
#
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
#
RUN pip install --upgrade pip
RUN pip install jupyterlab \
    pandas \
    numpy \
    pyyaml \
    Pillow \
    scipy \
    geopandas \
    descartes \
    rasterio \
    tensorflow-probability==0.7.0 \
    keras

RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" > \
	/etc/apt/sources.list.d/cran.list
# note the proxy for gpg
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
#
ENV DEBIAN_FRONTEND noninteractive
USER root

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
RUN apt-get update
#install R base dev
# Now install R and littler, and create a link for littler in /usr/local/bin
# Default CRAN repo is now set by R itself, and littler knows about it too
# r-cran-docopt is not currently in c2d4u so we install from source
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
	 littler \
	 r-base \
	 r-base-dev \
	 r-recommended \
	 python3-venv \
    && ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
    && ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
    && ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
    && ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
    && install.r docopt \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*

# install R packages
RUN R -e 'install.packages(c("tensorflow",\
                             "keras",\
                             "h2o",\
                             "tidyverse",\
                             "igraph",\
                             "DiagrammeR",\
                             "IRkernel",\
                             "abind",\
                             "moments",\
                             "mvtnorm",\
                             "ggthemes",\
                             "bayesplot",\
                             "HDInterval",\
                             "txtplot",\
                             "webshot"), dependencies=TRUE)'
RUN R -e 'IRkernel::installspec(user=FALSE)'
# here we use greta-dev instead of greta, because the library is still in development
# can use below once a tagged version is complete and on CRAN
# RUN R -e 'install.packages("greta")'
RUN Rscript -e "update.packages(ask = FALSE)" && \
    R -e "devtools::install_github('greta-dev/greta')"


#install RStudio
RUN printf '\nexport QT_STYLE_OVERRIDE=gtk\n' | tee -a ~/.profile

RUN update-locale
#old file location https://download2.rstudio.org/rstudio-server-1.1.456-amd64.deb
RUN wget  https://download2.rstudio.org/server/bionic/amd64/rstudio-server-1.2.1335-amd64.deb \
                                              && gdebi -n rstudio-server-1.2.1335-amd64.deb

RUN ln -s /usr/lib/rstudio-server/bin/pandoc/pandoc /usr/local/bin \
  && ln -s /usr/lib/rstudio-server/bin/pandoc/pandoc-citeproc /usr/local/bin \
  && git clone https://github.com/jgm/pandoc-templates \
  && mkdir -p /opt/pandoc/templates \
  && cp -r pandoc-templates*/* /opt/pandoc/templates && rm -rf pandoc-templates* \
  && mkdir /root/.pandoc && ln -s /opt/pandoc/templates /root/.pandoc/templates \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/ \
  && mkdir -p /usr/local/lib/R/etc \
  ## RStudio wants an /etc/R, will populate from $R_HOME/etc
  && mkdir -p /etc/R \
  ## Write config files in $R_HOME/etc
  && echo '\n\
    \n# Configure httr to perform out-of-band authentication if HTTR_LOCALHOST \
    \n# is not set since a redirect to localhost may not work depending upon \
    \n# where this Docker container is running. \
    \nif(is.na(Sys.getenv("HTTR_LOCALHOST", unset=NA))) { \
    \n  options(httr_oob_default = TRUE) \
    \n}' >> /usr/local/lib/R/etc/Rprofile.site \
  && echo "PATH=${PATH}" >> /usr/local/lib/R/etc/Renviron \
  ## Prevent rstudio from deciding to use /usr/bin/R if a user apt-get installs a package
  &&  echo 'rsession-which-r=/usr/local/bin/R' >> /etc/rstudio/rserver.conf \
  ## use more robust file locking to avoid errors when using shared volumes:
  && echo 'lock-type=advisory' >> /etc/rstudio/file-locks \
  ## configure git not to request password each time
  && git config --system credential.helper 'cache --timeout=3600' \
  && git config --system push.default simple \
  ## Set up S6 init system
  && wget -P /tmp/ https://github.com/just-containers/s6-overlay/releases/download/v1.11.0.1/s6-overlay-amd64.tar.gz \
  && tar xzf /tmp/s6-overlay-amd64.tar.gz -C / \
  && mkdir -p /etc/services.d/rstudio

ENV PATH=/usr/lib/rstudio-server/bin:$PATH

RUN echo '#!/usr/bin/with-contenv bash \
           \n exec /usr/lib/rstudio-server/bin/rserver --server-daemonize 0' \
           > /etc/services.d/rstudio/run \
   && echo '#!/bin/bash \
           \n rstudio-server stop' \
           > /etc/services.d/rstudio/finish

RUN wget https://github.com/rocker-org/rocker-versioned/blob/master/rstudio/userconf.sh
RUN cp userconf.sh /etc/cont-init.d/userconf

RUN wget https://github.com/rocker-org/rocker-versioned/blob/master/rstudio/pam-helper.sh
RUN cp pam-helper.sh /usr/lib/rstudio-server/bin/pam-helper

## running with "-e ADD=shiny" adds shiny server
#COPY add_shiny.sh /etc/cont-init.d/add

RUN ln /usr/bin/R /usr/local/bin/R

# Path of reticulate for Python (otherwise it looks for miniconda by default).
ENV RETICULATE_PYTHON="/usr/lib/python3.6"

EXPOSE 8787

## Need to configure non-root user for RStudio
RUN useradd rstudio \
    && echo "rstudio:rstudio" | chpasswd \
    && mkdir /home/rstudio \
    && chown rstudio:rstudio /home/rstudio \
    && addgroup rstudio staff

CMD ["/usr/lib/rstudio-server/bin/rserver", "--server-daemonize=0", "--server-app-armor-enabled=0"]
