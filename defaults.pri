CONFIG -= app_bundle
CONFIG -= qt

LIBS += -llapack -larmadillo -lhdf5 -lhdf5_cpp -lfftw3 -lfftw3_omp -lfftw3_threads
LIBS += -lyaml-cpp -lboost_system -lboost_filesystem
LIBS += -L/usr/local/lib -lopencv_core -lopencv_highgui
LIBS += -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_cpp

CONFIG += c++14
QMAKE_CXXFLAGS += -std=c++14  -fext-numeric-literals
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$QMAKE_CXXFLAGS -O3 -DARMA_NO_DEBUG

CURRENT_COMPILER = $$QMAKE_CXX
QMAKE_CXX = ccache $$CURRENT_COMPILER

INCLUDEPATH += /usr/local/include/opencv
INCLUDEPATH +=/usr/include/hdf5/serial/

INCLUDEPATH += $$TOP_PWD/src
SRC_DIR = $$TOP_PWD




