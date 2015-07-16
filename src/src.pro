include(../defaults.pri)

CONFIG   += console
CONFIG   -= app_bundle
CONFIG   -= qt
CONFIG += c++11

TEMPLATE = lib

TARGET = edog

SOURCES += \
    integrator.cpp \
    trapezoidal.cpp \
    outputmanager.cpp \
    lib.cpp \
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
    temporalKernels/diracDelta.cpp
HEADERS += \
    integrator.h \
    trapezoidal.h \
    outputmanager.h \
    lib.h \
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
    temporalKernels/diracDelta.h