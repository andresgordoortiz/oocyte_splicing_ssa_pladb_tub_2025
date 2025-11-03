# Use a base R image
FROM rocker/r-ver:4.3.1

# Install system dependencies for R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    pandoc \
    git \
    libxml2 \
    libxt6 \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre3-dev \
    libicu-dev \
    libjpeg-dev \
    libpng-dev \
    libglpk-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set environment variables
ENV RENV_VERSION=0.17.3
ENV R_LIBS_USER=/renv/library

# Copy renv.lock and other files to the container

COPY renv.lock /renv.lock
WORKDIR /
# Install renv
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')"


# Restore the R environment using renv
RUN R -e "tryCatch(renv::restore(), error = function(e) { Sys.sleep(10); renv::restore() })"

# By not setting an ENTRYPOINT, this Docker container is now ready to run any R script or RMarkdown file downstream.
