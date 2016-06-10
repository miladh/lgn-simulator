include(../../defaults.pri)
include(../apps_defaults.pri)

TARGET = lgnSimulator_stimuliAnalysis

SOURCES = stimulianalysis.cpp

OTHER_FILES += ./stimuliAnalysis.yaml

DISTFILES += \
    stimuliAnalysis.yaml
