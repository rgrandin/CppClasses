QT       += core

QT       -= gui

DEFINES += USEQT

TARGET = cpp-classes-test
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

include( ../src/include-cppclasses.pri )




HEADERS += \
    vtkxmliotest.h \
    xdmfiotest.h

SOURCES += main.cpp


CONFIG(debug, debug|release) {
    DESTDIR = debug
} else {
    DESTDIR = release
}



macx:{
    DEFINES += COMPILEMAC

    # OpenMP support -- GCC/MinGW
    #QMAKE_CXXFLAGS += -fopenmp
    #QMAKE_LFLAGS += -fopenmp

    # Enable C++11 support
    #QMAKE_CXXFLAGS += -std=c++11
    #DEFINES += CXX11

    #DEFINES += FFTW_TRANSPOSE
    #LIBS += -lfftw3 \       # Double-precision routines
    #        -lfftw3f        # Single-precision routines

}
unix:!macx{
    # "unix" is also valid for Mac system, so "!macx" is required for Linux-only

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

    # VTK libraries, version 5.10.1
    #INCLUDEPATH += /usr/local/include/vtk-5.10
    #LIBS += -L/usr/local/lib/vtk-5.10 -lvtkIO -lvtkCommon -lvtkImaging -lvtkFiltering -lvtkDICOMParser \
    #        -lvtkNetCDF -lvtkNetCDF_cxx -lvtkmetaio -lvtksqlite -lvtkpng -lvtkzlib \
    #        -lvtkjpeg -lvtkexpat -lvtksys -lvtkhdf5_hl -lvtkhdf5 -lvtktiff \
    #        -lLSDyna

    # VTK libraries, version 6.0
    INCLUDEPATH += /usr/local/include/vtk-6.0
    LIBS += -L/usr/local/lib -lvtkIOXML-6.0 -lvtkCommonCore-6.0 -lvtkIOImage-6.0


    # XDMF
    INCLUDEPATH += /usr/local/include                       # XDMF classes
    INCLUDEPATH += /usr/lib64/mpi/gcc/openmpi/include       # MPI header
    LIBS += -L/usr/local/xdmf/lib/ -lmetis -lvtkexoIIc -lvtkhdf5 -lvtklibxml2 -lXdmf -lXdmfUtils
    LIBS += -L/usr/lib64/mpi/gcc/openmpi/lib64 -lmpi -lmpi_cxx

}
win32:{
    DEFINES += COMPILEWIN32

    DEFINES += FFTW_TRANSPOSE
    INCLUDEPATH += C:\\FFTW
    QMAKE_LIBDIR += C:\\FFTW
    LIBS += -llibfftw3-3 -llibfftw3f-3 -llibfftw3l-3

    DEFINES += _CRT_SECURE_NO_WARNINGS
    DEFINES += CXX11
}

