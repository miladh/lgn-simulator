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
    fftIntegrator/integratortest_cosine.cpp \
    fftIntegrator/integratortest_patchgrating.cpp \
    fftIntegrator/integratortest_dampedsinusoid.cpp \
    fftIntegrator/integratortest_constant.cpp \
    fftIntegrator/integratortest_exponential.cpp

LIBS += -lunittest++ -L../src -ledog

DISTFILES += \
    configTests.cfg

HEADERS +=
