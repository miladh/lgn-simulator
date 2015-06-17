include(../defaults.pri)

CONFIG   += console
CONFIG   -= app_bundle
CONFIG   -= qt
CONFIG += c++11

TEMPLATE = lib

TARGET = edog

SOURCES += \
    impulseResponse.cpp \
    stimuli.cpp \
    integrator.cpp \
    trapezoidal.cpp \
    response.cpp \
    outputmanager.cpp \
    extendeddog.cpp \
    lib.cpp
HEADERS += \
    impulseResponse.h \
    stimuli.h \
    integrator.h \
    trapezoidal.h \
    response.h \
    outputmanager.h \
    extendeddog.h \
    lib.h
