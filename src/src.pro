include(../defaults.pri)

CONFIG   += console
CONFIG   -= app_bundle
CONFIG   -= qt
CONFIG += c++11

TEMPLATE = lib

TARGET = edog

SOURCES += \
    impulseResponse.cpp \
    integrator.cpp \
    trapezoidal.cpp \
    response.cpp \
    outputmanager.cpp \
    extendeddog.cpp \
    lib.cpp \
    ganglion/ganglion.cpp \
    ganglion/gangliondog.cpp \
    math/dog.cpp \
    stimuli/stimuli.cpp \
    stimuli/patchgrating.cpp \
    math/functions.cpp
HEADERS += \
    impulseResponse.h \
    integrator.h \
    trapezoidal.h \
    response.h \
    outputmanager.h \
    extendeddog.h \
    lib.h \
    ganglion/ganglion.h \
    ganglion/gangliondog.h \
    math/dog.h \
    stimuli/stimuli.h \
    stimuli/patchgrating.h \
    math/functions.h
