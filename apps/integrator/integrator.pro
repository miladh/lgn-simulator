include(../../defaults.pri)
include(../apps_defaults.pri)

TARGET = lgnSimulator_integrator

SOURCES += \
    integrator.cpp

OTHER_FILES += ./integrator.yaml

DISTFILES += \
    integrator.yaml
