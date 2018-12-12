FROM ubuntu:14.04

ENV HOME /root
WORKDIR /root

# Set the locale
# this is required so that ASPECT can format numbers with , separators
RUN locale-gen en_US.UTF-8  
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8

RUN apt-get update && apt-get -yq install gcc \
    g++ \
    gfortran \
    libblas-dev \
    liblapack-dev \
    libopenmpi-dev \
    openmpi-bin \
    cmake \
    git \
    wget \
    python

RUN wget http://www.netlib.org/scalapack/scalapack-2.0.2.tgz && tar -xzf scalapack-2.0.2.tgz
RUN cd scalapack-2.0.2 && cp SLmake.inc.example SLmake.inc && make
RUN cp scalapack-2.0.2/lib*.a /usr/local/lib/
COPY test_example .
RUN cd matrix_library && make && make test
