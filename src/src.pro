include(../defaults.pri)

CONFIG   += console
CONFIG   -= app_bundle
CONFIG   -= qt
CONFIG += c++11

TEMPLATE = lib

TARGET = edog


#LIBS += -llapack -larmadillo -lconfig++ -lboost_regex -lhdf5 -lhdf5_cpp -lfftw3 -lfftw3_omp -lfftw3_threads


#INCLUDEPATH += /usr/include/hdf5/serial

SOURCES += \
    integrator.cpp \
    outputmanager.cpp \
    stimuli/stimuli.cpp \
    math/functions.cpp \
    neurons/neuron.cpp \
    neurons/corticalcell.cpp \
    temporalKernels/temporalkernel.cpp \
    spatialKernels/spatialkernel.cpp \
    neurons/relaycell.cpp \
    spatialKernels/dog.cpp \
    temporalKernels/decayingexponential.cpp \
    neurons/ganglioncell.cpp \
    spatialKernels/ellipticgaussian.cpp \
    temporalKernels/dampedoscillator.cpp \
    math/ffthelper.cpp \
    neurons/interneuron.cpp \
    temporalKernels/temporallyconstant.cpp \
    stimuli/naturalscene.cpp \
    stimuli/staticimage.cpp \
    stimuli/naturalscenevideo.cpp \
    temporalKernels/temporaldelta.cpp \
    stimuli/grating/grating.cpp \
    stimuli/grating/gaussianmaskgrating.cpp \
    stimuli/grating/circlemaskgrating.cpp \
    stimuli/grating/fullfieldgrating.cpp \
    staticNonlinearity/staticnonlinearity.cpp \
    staticNonlinearity/thresholdnonlinearity.cpp \
    staticNonlinearity/heavisidenonlinearity.cpp \
    staticNonlinearity/sigmoidalnonlinearity.cpp
HEADERS += \
    integrator.h \
    outputmanager.h \
    stimuli/stimuli.h \
    math/functions.h \
    neurons/neuron.h \
    neurons/corticalcell.h \
    temporalKernels/temporalkernel.h \
    spatialKernels/spatialkernel.h \
    neurons/relaycell.h \
    spatialKernels/dog.h \
    temporalKernels/decayingexponential.h \
    neurons/ganglioncell.h \
    spatialKernels/ellipticgaussian.h \
    temporalKernels/dampedoscillator.h \
    math/ffthelper.h \
    neurons/interneuron.h \
    temporalKernels/temporallyconstant.h \
    stimuli/naturalscene.h \
    stimuli/staticimage.h \
    stimuli/naturalscenevideo.h \
    temporalKernels/temporaldelta.h \
    stimuli/grating/grating.h \
    stimuli/grating/gaussianmaskgrating.h \
    stimuli/grating/circlemaskgrating.h \
    stimuli/grating/fullfieldgrating.h \
    staticNonlinearity/staticnonlinearity.h \
    staticNonlinearity/thresholdnonlinearity.h \
    staticNonlinearity/heavisidenonlinearity.h \
    staticNonlinearity/sigmoidalnonlinearity.h
