include(../defaults.pri)

TEMPLATE = lib
TARGET = ../lib/edog


SOURCES += \
    integrator.cpp \
    outputmanager.cpp \
    stimuli/stimuli.cpp \
    math/functions.cpp \
    neurons/neuron.cpp \
    neurons/corticalcell.cpp \
    neurons/relaycell.cpp \
    neurons/ganglioncell.cpp \
    math/ffthelper.cpp \
    neurons/interneuron.cpp \
    stimuli/naturalscene.cpp \
    stimuli/staticimage.cpp \
    stimuli/naturalscenevideo.cpp \
    stimuli/grating/grating.cpp \
    stimuli/grating/gaussianmaskgrating.cpp \
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
    kernels/spatialKernels/gaussian.cpp \
    kernels/spatialKernels/spatialdelta.cpp \
    kernels/spatialKernels/spatialkernel.cpp \
    kernels/spatialKernels/spatiallyconstant.cpp \
    kernels/temporalKernels/dampedoscillator.cpp \
    kernels/temporalKernels/decayingexponential.cpp \
    kernels/temporalKernels/temporaldelta.cpp \
    kernels/temporalKernels/temporalGaussian.cpp \
    kernels/temporalKernels/temporalkernel.cpp \
    kernels/temporalKernels/temporallyconstant.cpp

HEADERS += \
    integrator.h \
    outputmanager.h \
    stimuli/stimuli.h \
    math/functions.h \
    neurons/neuron.h \
    neurons/corticalcell.h \
    neurons/relaycell.h \
    neurons/ganglioncell.h \
    math/ffthelper.h \
    neurons/interneuron.h \
    stimuli/naturalscene.h \
    stimuli/staticimage.h \
    stimuli/naturalscenevideo.h \
    stimuli/grating/grating.h \
    stimuli/grating/gaussianmaskgrating.h \
    stimuli/grating/circlemaskgrating.h \
    stimuli/grating/fullfieldgrating.h \
    staticNonlinearity/staticnonlinearity.h \
    staticNonlinearity/thresholdnonlinearity.h \
    staticNonlinearity/heavisidenonlinearity.h \
    staticNonlinearity/sigmoidalnonlinearity.h \
    ../include/lgnSimulator.h \
    kernels/kernel.h \
    kernels/separablekernel.h \
    kernels/temporalKernels/dampedoscillator.h \
    kernels/temporalKernels/decayingexponential.h \
    kernels/temporalKernels/temporaldelta.h \
    kernels/temporalKernels/temporalGaussian.h \
    kernels/temporalKernels/temporalkernel.h \
    kernels/temporalKernels/temporallyconstant.h \
    kernels/spatialKernels/dog.h \
    kernels/spatialKernels/ellipticgaussian.h \
    kernels/spatialKernels/gaussian.h \
    kernels/spatialKernels/spatialdelta.h \
    kernels/spatialKernels/spatialkernel.h \
    kernels/spatialKernels/spatiallyconstant.h

OTHER_FILES +=

!equals(PWD, $${OUT_PWD}) {
    QMAKE_POST_LINK += $(COPY_DIR) $$OUT_PWD/../lib $$PWD/../
}
