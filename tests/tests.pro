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
    integratorTests/test_integrator_constant.cpp \
    integratorTests/test_integrator_grating.cpp \
    fftHelperTests/test_helper_fftfreq.cpp \
    fftHelperTests/test_helper_fftshift1d.cpp \
    fftHelperTests/test_helper_fftshift2d.cpp \
    fftHelperTests/test_helper_fftshift3d.cpp \
    fftHelperTests/test_helper_ifftshift3d.cpp \
    integratorTests/convolutionTests/test_conv_gaussdelta.cpp \
    integratorTests/convolutionTests/test_conv_gaussconstant.cpp \
    integratorTests/convolutionTests/test_conv_dampedoscdelta.cpp \
    integratorTests/convolutionTests/test_conv_decayingExpDelta.cpp




HEADERS += \
    systemTests/kernelsettings.h
