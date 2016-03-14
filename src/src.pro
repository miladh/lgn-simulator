include(../defaults.pri)

TEMPLATE = lib
TARGET = ../lib/lgn-simulator


SOURCES += \
    integrator.cpp \
    outputmanager.cpp \
    neurons/neuron.cpp \
    neurons/corticalcell.cpp \
    neurons/relaycell.cpp \
    neurons/ganglioncell.cpp \
    neurons/interneuron.cpp \
    stimuli/naturalscene.cpp \
    stimuli/staticimage.cpp \
    stimuli/naturalscenevideo.cpp \
    stimuli/grating/grating.cpp \
    stimuli/grating/circlemaskgrating.cpp \
    stimuli/grating/fullfieldgrating.cpp \
    staticNonlinearity/staticnonlinearity.cpp \
    staticNonlinearity/thresholdnonlinearity.cpp \
    staticNonlinearity/heavisidenonlinearity.cpp \
    staticNonlinearity/sigmoidalnonlinearity.cpp \
    kernels/kernel.cpp \
    kernels/separablekernel.cpp \
    kernels/spatialKernels/dog.cpp \
    kernels/spatialKernels/ellipticgaussian.cpp \
    kernels/spatialKernels/spatialdelta.cpp \
    kernels/spatialKernels/spatialkernel.cpp \
    kernels/spatialKernels/spatiallyconstant.cpp \
    kernels/temporalKernels/decayingexponential.cpp \
    kernels/temporalKernels/temporaldelta.cpp \
    kernels/temporalKernels/temporalGaussian.cpp \
    kernels/temporalKernels/temporalkernel.cpp \
    kernels/temporalKernels/temporallyconstant.cpp \
    helper/ffthelper.cpp \
    helper/special.cpp \
    stimuli/stimulus.cpp \
    kernels/temporalKernels/twosidedexponentialdecay.cpp \
    kernels/temporalKernels/doe.cpp \
    kernels/temporalKernels/biphasic.cpp \
    kernels/spatialKernels/spatialgaussian.cpp

HEADERS += \
    integrator.h \
    outputmanager.h \
    neurons/neuron.h \
    neurons/corticalcell.h \
    neurons/relaycell.h \
    neurons/ganglioncell.h \
    neurons/interneuron.h \
    stimuli/naturalscene.h \
    stimuli/staticimage.h \
    stimuli/naturalscenevideo.h \
    stimuli/grating/grating.h \
    stimuli/grating/circlemaskgrating.h \
    stimuli/grating/fullfieldgrating.h \
    staticNonlinearity/staticnonlinearity.h \
    staticNonlinearity/thresholdnonlinearity.h \
    staticNonlinearity/heavisidenonlinearity.h \
    staticNonlinearity/sigmoidalnonlinearity.h \
    ../include/lgnSimulator.h \
    kernels/kernel.h \
    kernels/separablekernel.h \
    kernels/temporalKernels/decayingexponential.h \
    kernels/temporalKernels/temporaldelta.h \
    kernels/temporalKernels/temporalGaussian.h \
    kernels/temporalKernels/temporalkernel.h \
    kernels/temporalKernels/temporallyconstant.h \
    kernels/spatialKernels/dog.h \
    kernels/spatialKernels/ellipticgaussian.h \
    kernels/spatialKernels/spatialdelta.h \
    kernels/spatialKernels/spatialkernel.h \
    kernels/spatialKernels/spatiallyconstant.h \
    helper/ffthelper.h \
    helper/helperconstants.h \
    helper/special.h \
    stimuli/stimulus.h \
    kernels/temporalKernels/twosidedexponentialdecay.h \
    kernels/temporalKernels/doe.h \
    kernels/temporalKernels/biphasic.h \
    kernels/spatialKernels/spatialgaussian.h

OTHER_FILES +=

!equals(PWD, $${OUT_PWD}) {
    QMAKE_POST_LINK += $(COPY_DIR) $$OUT_PWD/../lib $$PWD/../
}
