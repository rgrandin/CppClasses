INCLUDEPATH += $$PWD

HEADERS += \
    $$PWD/ProgressInfo.h \
    $$PWD/QtIntermediary.h \
    $$PWD/QtIntermediaryBase.h \
    $$PWD/ArrayBase.h \
    $$PWD/Array1D.h \
    $$PWD/Array2D.h \
    $$PWD/Array3D.h \
    $$PWD/Array4D.h \
    $$PWD/PArrayBase.h \
    $$PWD/PArray1D.h \
    $$PWD/PArray2D.h \
    $$PWD/PArray3D.h \
    $$PWD/PArray4D.h \
    $$PWD/DataFilters.h \
    $$PWD/DLList.h \
    $$PWD/SLList.h \
    $$PWD/QSafeApplication.h \
    $$PWD/Stats.h \
    $$PWD/StringManip.h \
    $$PWD/UniformVolume.h \
    $$PWD/UnitConvert.h \
    $$PWD/XdmfIO.h

SOURCES += \
    $$PWD/ProgressInfo.cpp \
    $$PWD/QtIntermediary.cpp \
    $$PWD/QtIntermediaryBase.cpp \
    $$PWD/StringManip.cpp



# ===============================================================================
#
# Add source files for templated classes.  Immediately remove them from 'SOURCES'
# so that they are not compiled (templates must be all in the headers so the
# compiler can create an appropriate instance at compile-time).
#
#

SOURCES += \
    $$PWD/ArrayBase.cpp \
    $$PWD/Array1D.cpp \
    $$PWD/Array2D.cpp \
    $$PWD/Array3D.cpp \
    $$PWD/Array4D.cpp \
    $$PWD/PArrayBase.cpp \
    $$PWD/PArray1D.cpp \
    $$PWD/PArray2D.cpp \
    $$PWD/PArray3D.cpp \
    $$PWD/PArray4D.cpp \
    $$PWD/DataFilters.cpp \
    $$PWD/DLList.cpp \
    $$PWD/SLList.cpp \
    $$PWD/Stats.cpp \
    $$PWD/UniformVolume.cpp \
    $$PWD/UnitConvert.cpp

SOURCES -= \
    $$PWD/ArrayBase.cpp \
    $$PWD/Array1D.cpp \
    $$PWD/Array2D.cpp \
    $$PWD/Array3D.cpp \
    $$PWD/Array4D.cpp \
    $$PWD/PArrayBase.cpp \
    $$PWD/PArray1D.cpp \
    $$PWD/PArray2D.cpp \
    $$PWD/PArray3D.cpp \
    $$PWD/PArray4D.cpp \
    $$PWD/DataFilters.cpp \
    $$PWD/DLList.cpp \
    $$PWD/SLList.cpp \
    $$PWD/Stats.cpp \
    $$PWD/UniformVolume.cpp \
    $$PWD/UnitConvert.cpp
