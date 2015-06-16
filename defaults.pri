INCLUDEPATH += $$PWD/src
SRC_DIR = $$PWD

LIBS += -llapack -larmadillo -lconfig++ -lhdf5 -lhdf5_cpp

QMAKE_CXXFLAGS += -std=c++0x
