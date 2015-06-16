TEMPLATE = subdirs
CONFIG+=ordered


SUBDIRS = \
    src \
    app \
    tests

app.depends = src
tests.depends = src

OTHER_FILES += \
    defaults.pri

LIBS += -llapack -larmadillo -lconfig++ -lboost_regex -lhdf5 -lhdf5_cpp

