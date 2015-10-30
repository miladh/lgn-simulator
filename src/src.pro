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
    spatialKernels/ellipticgaussian.cpp \
    temporalKernels/dampedoscillator.cpp \
    math/ffthelper.cpp \
    stimuli/grating.cpp \
    neurons/interneuron.cpp \
    stimuli/oscillatinggaussian.cpp \
    temporalKernels/temporallyconstant.cpp \
    stimuli/naturalscene.cpp \
    stimuli/staticimage.cpp
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
    spatialKernels/ellipticgaussian.h \
    temporalKernels/dampedoscillator.h \
    math/ffthelper.h \
    stimuli/grating.h \
    neurons/interneuron.h \
    stimuli/oscillatinggaussian.h \
    temporalKernels/temporallyconstant.h \
    stimuli/naturalscene.h \
    stimuli/staticimage.h
