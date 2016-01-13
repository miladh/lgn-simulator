INCLUDEPATH += $$PWD/src
SRC_DIR = $$PWD
#INCLUDE_DIRS += /usr/include/hdf5/serial/
LIBS += -llapack -larmadillo -lconfig++ -lhdf5 -lhdf5_cpp

INCLUDEPATH += /usr/local/include/opencv
INCLUDEPATH +=/usr/include/hdf5/serial/



LIBS += -L/usr/local/lib -lopencv_core -lopencv_highgui
LIBS += -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_cpp

QMAKE_CXXFLAGS += -std=c++14

CONFIG += c++14


LIBS += -llapack -larmadillo -lboost_regex -lhdf5 -lhdf5_cpp -lfftw3 -lfftw3_omp -lfftw3_threads
LIBS += -lyaml-cpp
