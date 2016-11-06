FROM ubuntu:16.04

RUN apt-get -qq update
RUN apt-get install -qq sudo

RUN sudo apt-get -qq update
RUN sudo apt-get -qq install software-properties-common
RUN sudo apt-get -qq install python-software-properties
RUN sudo apt-add-repository -y ppa:beineri/opt-qt57-xenial
RUN sudo apt-get -qq update
RUN sudo apt-get -qq install qt573d qt57charts-no-lgpl qt57datavis-no-lgpl qt57declarative qt57svg


# needed by build
RUN sudo apt-get -qq install libfftw3-dev\
    libhdf5-dev\
    libyaml-dev\
    libyaml-cpp-dev\
    libgsl-dev\
    catch\
    libboost-dev\
    libarmadillo-dev\
    python-yaml


WORKDIR /app

ADD . /app/
