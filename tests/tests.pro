include(../defaults.pri)

CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TEMPLATE = app

SOURCES += main.cpp \
    integratortests.cpp \
    stimulitests.cpp \
    mathTests.cpp \
    fft1DTests.cpp \
    fft2DTests.cpp \
    ffthelpertests.cpp

LIBS += -lunittest++ -L../src -ledog

DISTFILES += \
    configTests.cfg
