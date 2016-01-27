include(../../defaults.pri)
include(../apps_defaults.pri)

TARGET = edog_nonlaggedXCells

SOURCES = nonlaggedXCells.cpp

OTHER_FILES += ./nonlaggedXCellsConfig.yaml

DISTFILES += \
    nonlaggedXCellsConfig.yaml
