include(../defaults.pri)

CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TEMPLATE = app

SOURCES += main.cpp \
    integratortests.cpp \
    fft.cpp

LIBS += -lunittest++ -L../src -ledog

DISTFILES += \
    configTests.cfg
