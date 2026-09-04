FROM docker.io/library/fedora@sha256:be9d65e2344d805cc11114319c685ecaa96b6d9b4350a0a6460cdb931babbd19

RUN dnf -y install git subversion make gcc gcc-gfortran gcc-c++ tar xz curl which \
    libcurl-devel openssl-devel zlib-devel bzip2-devel xz-devel pcre2-devel readline-devel \
    libicu-devel libjpeg-turbo-devel libpng-devel libtiff-devel cairo-devel texinfo perl \
 && git clone --depth 1 --branch 2.0 https://github.com/flame/blis.git /opt/blis \
 && cd /opt/blis && git fetch --depth 1 origin e8566eb3e773fb54d11b33e371d13f22d2941e50 \
 && git checkout --detach e8566eb3e773fb54d11b33e371d13f22d2941e50 \
 && ./configure auto --enable-threading=no --prefix=/opt/blis-install \
 && make -j1 && make install \
 && svn export -r 90483 https://svn.r-project.org/R/trunk /opt/R-src \
 && cd /opt/R-src && ./configure --with-blas='-L/opt/blis-install/lib -lblis' --prefix=/opt/R/devel-blis \
 && make -j1 && make install \
 && /opt/R/devel-blis/bin/R --vanilla -s -e 'dir.create("/opt/blis-library", recursive=TRUE); install.packages(c("testthat", "Rfast", "matrixStats", "knitr", "rmarkdown", "ashr", "MASS", "susieR"), lib="/opt/blis-library", repos="https://cloud.r-project.org", dependencies=c("Depends", "Imports", "LinkingTo", "Suggests"), type="source")' \
 && printf '%s\n' '#!/bin/sh' 'set -eu' 'exec /opt/R/devel-blis/bin/R CMD check --no-manual --no-build-vignettes "$1/"*.tar.gz' > /usr/local/bin/r-check \
 && chmod 0755 /usr/local/bin/r-check \
 && dnf clean all

ENV BLIS_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 R_LIBS_USER=/opt/blis-library PATH=/opt/R/devel-blis/bin:$PATH LD_LIBRARY_PATH=/opt/blis-install/lib
