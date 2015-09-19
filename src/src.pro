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
    stimuli/patchgrating.cpp \
    math/functions.cpp \
    neurons/neuron.cpp \
    neurons/corticalcell.cpp \
    temporalKernels/temporalkernel.cpp \
    spatialKernels/spatialkernel.cpp \
    neurons/relaycell.cpp \
    spatialKernels/dog.cpp \
    temporalKernels/decayingexponential.cpp \
    spatialKernels/gaussian.cpp \
    neurons/ganglioncell.cpp \
    temporalKernels/diracDelta.cpp \
    spatialKernels/ellipticgaussian.cpp \
    temporalKernels/dampedoscillator.cpp \
    stimuli/dogstim.cpp \
    integratorsettings.cpp
HEADERS += \
    integrator.h \
    outputmanager.h \
    stimuli/stimuli.h \
    stimuli/patchgrating.h \
    math/functions.h \
    neurons/neuron.h \
    neurons/corticalcell.h \
    temporalKernels/temporalkernel.h \
    spatialKernels/spatialkernel.h \
    neurons/relaycell.h \
    spatialKernels/dog.h \
    temporalKernels/decayingexponential.h \
    spatialKernels/gaussian.h \
    neurons/ganglioncell.h \
    temporalKernels/diracDelta.h \
    spatialKernels/ellipticgaussian.h \
    temporalKernels/dampedoscillator.h \
    stimuli/dogstim.h \
    integratorsettings.h
