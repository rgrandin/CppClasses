QT       += core

QT       -= gui

DEFINES += USEQT

TARGET = cpp-classes-test
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

include( ../src/include-cppclasses.pri )

#HEADERS +=  \
#            ../src/QtIntermediaryBase.h

#SOURCES += main.cpp \
#            ../src/QtIntermediaryBase.cpp

SOURCES += main.cpp


CONFIG(debug, debug|release) {
    DESTDIR = debug
} else {
    DESTDIR = release
}




unix:{
    DEFINES += COMPILELINUX

    # OpenMP support -- GCC/MinGW
    QMAKE_CXXFLAGS += -fopenmp
    QMAKE_LFLAGS += -fopenmp

    # Enable C++11 support
    QMAKE_CXXFLAGS += -std=c++11
    DEFINES += CXX11

    DEFINES += FFTW_TRANSPOSE
    LIBS += -lfftw3 \       # Double-precision routines
            -lfftw3f        # Single-precision routines

}
win32:{
    DEFINES += FFTW_TRANSPOSE
    INCLUDEPATH += C:\\FFTW
    QMAKE_LIBDIR += C:\\FFTW
    LIBS += -llibfftw3-3 -llibfftw3f-3 -llibfftw3l-3

    DEFINES += _CRT_SECURE_NO_WARNINGS
    DEFINES += CXX11
}
