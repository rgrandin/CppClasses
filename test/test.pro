QT       += core

QT       -= gui

TARGET = cpp-classes-test
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp


CONFIG(debug, debug|release) {
    DESTDIR = debug
} else {
    DESTDIR = release
}


INCLUDEPATH += ../src


unix:{
    DEFINES += COMPILELINUX

    # OpenMP support -- GCC/MinGW
    QMAKE_CXXFLAGS += -fopenmp
    QMAKE_LFLAGS += -fopenmp

    # Enable C++11 support
    QMAKE_CXXFLAGS += -std=c++11
    DEFINES += CXX11
}
win32:{
    DEFINES += _CRT_SECURE_NO_WARNINGS
    DEFINES += CXX11
}
