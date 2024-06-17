FROM rocker/rstudio:4.3.1

RUN apt-get update && \
  apt-get install -y \
  libpng-dev \
  libcurl4-openssl-dev \
  libxml2-dev \
  libssl-dev \
  pandoc \
  libfontconfig1-dev \
  nano \
  libhdf5-dev \
  default-jdk \
  cmake \
  libfftw3-dev \
  libgeos-dev \
  libmagick++-dev \
  libproj-dev \
  libgdal-dev \
  libharfbuzz-dev \
  libfribidi-dev  \
  libudunits2-dev \
  libgsl-dev \
  libgmp3-dev \
  libglpk40 \
  patch \
  libglpk-dev \
  libgmp-dev

RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

WORKDIR /project
COPY renv.lock renv.lock

# approach one
ENV RENV_PATHS_LIBRARY renv/library

RUN R -e "renv::restore()"

# RUN R -e "renv::restore(repos = 'https://cloud.r-project.org/')"

# RUN R -e "options(repos = c(CRAN = 'https://cloud.r-project.org'))"
