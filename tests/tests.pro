include(../defaults.pri)

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TARGET = lgnSimulator_tests

LIBS += -L$$TOP_OUT_PWD/lib -llgn-simulator
INCLUDEPATH += $$TOP_PWD/include


SOURCES += main.cpp \
    test_special.cpp \
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
    integratorTests/test_integrator_patchgrating.cpp \
    integratorTests/test_integrator_gaussiangrating.cpp \
    systemTests/test_system_g_grating.cpp \
    systemTests/test_system_gi_grating.cpp \
    systemTests/test_system_gr_grating.cpp \
    systemTests/test_system_grc_grating.cpp \
    systemTests/test_system_gric_grating.cpp \
    stimulusTests/test_stimulus_circlegrating.cpp \
    stimulusTests/test_stimulus_fullfieldgrating.cpp \
    systemTests/test_system_g_patchgrating.cpp





HEADERS +=
