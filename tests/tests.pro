include(../defaults.pri)

CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TEMPLATE = app

SOURCES += main.cpp \
    stimulitests.cpp \
    mathTests.cpp \
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
    fftIntegrator/integratortest_constant.cpp \
    fftIntegrator/integratortest_exponential.cpp \
    systemTests/systemtest_gr.cpp \
    systemTests/systemtest_g.cpp \
    systemTests/kernelsettings.cpp \
    systemTests/systemtest_grc.cpp \
    systemTests/systemtest_gric.cpp \
    systemTests/systemtest_gi.cpp

LIBS += -lunittest++ -L../src -ledog

DISTFILES += \
    configTests.cfg

HEADERS += \
    systemTests/kernelsettings.h
