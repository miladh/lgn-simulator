include(../defaults.pri)

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TARGET = lgnSimulator_tests

LIBS += -lunittest++
LIBS += -L$$TOP_OUT_PWD/lib -llgn-simulator
INCLUDEPATH += $$TOP_PWD/include


SOURCES += main.cpp \
    fftHelper/fftfreqtests.cpp \
    fftHelper/fftshift1dtests.cpp \
    fftHelper/fftshift2dtests.cpp \
    fftHelper/fftshift3dtests.cpp \
    fftHelper/ifftshift3dtests.cpp \
    systemTests/systemtest_gr.cpp \
    systemTests/systemtest_g.cpp \
    systemTests/kernelsettings.cpp \
    systemTests/systemtest_grc.cpp \
    systemTests/systemtest_gric.cpp \
    systemTests/systemtest_gi.cpp \
    systemTests/systemtest_spot.cpp \
    test_special.cpp \
    stimulusTests/test_stimulus_fullfieldgrating.cpp \
    test_kernel.cpp \
    integratorTests/integratortest_exponential.cpp \
    integratorTests/test_integrator_constant.cpp \
    integratorTests/test_integrator_convolution.cpp \
    integratorTests/test_integrator_grating.cpp




HEADERS += \
    systemTests/kernelsettings.h
