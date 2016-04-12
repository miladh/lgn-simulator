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
    systemTests/kernelsettings.cpp \
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
    integratorTests/convolutionTests/test_conv_decayingexpdelta.cpp \
    integratorTests/convolutionTests/test_conv_deltadelta.cpp \
    integratorTests/convolutionTests/test_conv_biphasicdelta.cpp \
    integratorTests/convolutionTests/test_conv_biphasicdoe.cpp \
    systemTests/test_system_g.cpp \
    systemTests/test_system_gi.cpp \
    systemTests/test_system_gr.cpp \
    systemTests/test_system_grc.cpp \
    systemTests/test_system_gric.cpp \
    systemTests/test_system_g_spot.cpp




HEADERS += \
    systemTests/kernelsettings.h
