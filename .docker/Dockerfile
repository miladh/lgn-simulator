FROM ubuntu:16.04

RUN apt-get -qq update
RUN apt-get install -qq sudo

RUN sudo apt-get -qq update
RUN sudo apt-get -qq install software-properties-common
RUN sudo apt-get -qq install python-software-properties
RUN sudo apt-get -qq install qt5-default


# install packages needed by configure
RUN sudo apt-get -qq install patch git build-essential curl

# needed by build
RUN sudo apt-get -qq install ccache\
    libfftw3-dev\
    libhdf5-dev\
    libyaml-dev\
    libyaml-cpp-dev\
    libgsl-dev\
    catch\
    libboost-dev\
    libboost-system-dev\
    libboost-filesystem-dev\
    libarmadillo-dev\
    libopencv-dev
    
