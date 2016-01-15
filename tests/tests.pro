include(../defaults.pri)

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TARGET = edog-tests

LIBS += -lunittest++ -L$$TOP_OUT_PWD/lib -ledog
INCLUDEPATH += $$TOP_PWD/include


SOURCES += main.cpp \
    stimulitests.cpp \
    mathTests.cpp \
    fftHelper/fftfreqtests.cpp \
    fftHelper/fftshift1dtests.cpp \
    fftHelper/fftshift2dtests.cpp \
    fftHelper/fftshift3dtests.cpp \
    fftHelper/ifftshift3dtests.cpp \
    fftIntegrator/integratortest_dog.cpp \
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




HEADERS += \
    systemTests/kernelsettings.h
