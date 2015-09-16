INCLUDEPATH += $$PWD/src
SRC_DIR = $$PWD

LIBS += -llapack -larmadillo -lconfig++ -lhdf5 -lhdf5_cpp
INCLUDEPATH += /usr/local/include/opencv
LIBS += -L/usr/local/lib -lopencv_core -lopencv_highgui
QMAKE_CXXFLAGS += -std=c++0x
