FROM docker.io/library/fedora@sha256:be9d65e2344d805cc11114319c685ecaa96b6d9b4350a0a6460cdb931babbd19

RUN dnf -y install subversion make gcc gcc-gfortran gcc-c++ tar xz curl which \
    libcurl-devel openssl-devel zlib-devel bzip2-devel xz-devel pcre2-devel readline-devel \
    libicu-devel libjpeg-turbo-devel libpng-devel libtiff-devel cairo-devel texinfo perl \
 && svn export -r 90483 https://svn.r-project.org/R/trunk /opt/R-src \
 && cd /opt/R-src && ./configure --disable-openmp --prefix=/opt/R/devel-noomp \
 && make -j1 && make install \
 && /opt/R/devel-noomp/bin/R --vanilla -s -e 'dir.create("/opt/noomp-library", recursive=TRUE); install.packages(c("testthat", "Rfast", "matrixStats", "knitr", "rmarkdown", "ashr", "MASS", "susieR"), lib="/opt/noomp-library", repos="https://cloud.r-project.org", dependencies=c("Depends", "Imports", "LinkingTo", "Suggests"), type="source")' \
 && printf '%s\n' '#!/bin/sh' 'set -eu' 'exec /opt/R/devel-noomp/bin/R CMD check --no-manual --no-build-vignettes "$1/"*.tar.gz' > /usr/local/bin/r-check \
 && chmod 0755 /usr/local/bin/r-check \
 && dnf clean all

ENV SHLIB_OPENMP_CFLAGS= SHLIB_OPENMP_CXXFLAGS= SHLIB_OPENMP_FFLAGS= SHLIB_OPENMP_LDFLAGS= R_LIBS_USER=/opt/noomp-library PATH=/opt/R/devel-noomp/bin:$PATH
