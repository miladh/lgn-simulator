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
    fftHelper/fftfreqtests.cpp \
    fftHelper/fftshift1dtests.cpp \
    fftHelper/fftshift2dtests.cpp \
    fftHelper/fftshift3dtests.cpp \
    fftHelper/ifftshift3dtests.cpp

LIBS += -lunittest++ -L../src -ledog

DISTFILES += \
    configTests.cfg
