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
    ganglion/ganglion.cpp \
    ganglion/gangliondog.cpp \
    math/dog.cpp \
    stimuli/stimuli.cpp \
    stimuli/patchgrating.cpp \
    math/functions.cpp \
    relay/relay.cpp \
    relay/originaledog.cpp \
    neurons/neuron.cpp \
    neurons/corticalcell.cpp \
    temporalKernels/temporalkernel.cpp \
    spatialKernels/spatialkernel.cpp \
    neurons/relaycell.cpp
HEADERS += \
    integrator.h \
    trapezoidal.h \
    outputmanager.h \
    lib.h \
    ganglion/ganglion.h \
    ganglion/gangliondog.h \
    math/dog.h \
    stimuli/stimuli.h \
    stimuli/patchgrating.h \
    math/functions.h \
    relay/relay.h \
    relay/originaledog.h \
    neurons/neuron.h \
    neurons/corticalcell.h \
    temporalKernels/temporalkernel.h \
    spatialKernels/spatialkernel.h \
    neurons/relaycell.h
