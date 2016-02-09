include(../../defaults.pri)
include(../apps_defaults.pri)

TARGET = lgnSimulator_nonlaggedXCells

SOURCES = nonlaggedXCells.cpp

OTHER_FILES += ./nonlaggedXCellsConfig.yaml

DISTFILES += \
    nonlaggedXCellsConfig.yaml


