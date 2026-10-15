# rms_gene_explorer — runtime image.
#
# Base: rocker/shiny-verse pins both R and the CRAN snapshot to the R release date, so
# builds are reproducible. 4.5.3 matches the R the app was developed and tested against.
#
# shiny-verse already ships: shiny, bslib (a shiny dependency), ggplot2, dplyr, DBI,
# data.table, fst and arrow. Four packages the app needs are NOT in it, plus rprojroot,
# which arrives only indirectly via devtools/testthat and is too cheap to rely on.
#
# NOTE: every rocker/shiny* tag is published for amd64 only. On Apple Silicon this image
# runs under emulation — build and run it with `--platform linux/amd64`.
FROM rocker/shiny-verse:4.5.3

# ggiraph pulls systemfonts and gdtools, which need fontconfig and freetype present at
# runtime, not just at build time.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libfontconfig1-dev \
        libfreetype6-dev \
        libpng-dev \
        libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# --skipinstalled so anything shiny-verse already has is left alone.
RUN install2.r --error --skipinstalled --ncpus -1 \
        ggiraph \
        DT \
        duckdb \
        config \
        rprojroot

# bslib compiles its Bootstrap theme through sass and caches the result in
# <appdir>/app_cache. Shiny fixes that location, so it cannot be redirected by an env var —
# the run command mounts the code directory writable and the store read-only instead.

# The app reads only `store/`; the ETL scripts and the upstream source_dir are not needed
# at runtime, so nothing controlled-access has to be mounted beyond the store itself.
WORKDIR /app

EXPOSE 3838

# Run the app, NOT Shiny Server.
#
# Without this CMD the image inherits `CMD ["/init"]` from rocker/shiny-verse: the s6
# supervisor starts the shiny-server daemon, which serves /srv/shiny-server (rocker copies
# the shiny package's example apps there). A bare `docker run <image>` then shows
# "Welcome to Shiny Server!" and the bind-mounted repo is never read — the app looks
# broken when in fact it was never asked to start.
#
# The app path is RELATIVE so it resolves against whatever working directory the caller
# sets: `-w /app` with the repo mounted at /app, or `-w "$PWD"` with the repo mounted at
# its own host path. Both work.
#
# host='0.0.0.0' is the container's own interface and is required for the published port
# to reach it. Access control belongs on the host side of `-p` (127.0.0.1), not here.
CMD ["Rscript", "-e", "shiny::runApp('shiny/rms_gene_explorer', host = '0.0.0.0', port = 3838, launch.browser = FALSE)"]
