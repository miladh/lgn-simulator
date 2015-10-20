include(../defaults.pri)

CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TEMPLATE = app

SOURCES += main.cpp \
    stimulitests.cpp \
    mathTests.cpp \
    fft1DTests.cpp \
    fftHelper/fftfreqtests.cpp \
    fftHelper/fftshift1dtests.cpp \
    fftHelper/fftshift2dtests.cpp \
    fftHelper/fftshift3dtests.cpp \
    fftHelper/ifftshift3dtests.cpp \
    fftIntegrator/integratortest_dog.cpp \
    fftIntegrator/integratortest_gauss.cpp \
    fftIntegrator/integratortest_grating.cpp \
    fftIntegrator/integratortest_cosine.cpp

LIBS += -lunittest++ -L../src -ledog

DISTFILES += \
    configTests.cfg
