
# ==============================================================================
#
#   NOTES
#
# ==============================================================================
#
#	1.) These files are intended to be compiled into other projects and
#	    not as a stand-alone library.  Attempting to compile as a stand-alone
#	    library will cause errors due to "PARENT_SCOPE" variable options.
#	
#	    To incorporate this code into a project, copy the 'src' directory
#	    to the source directory of the project, and configure the project
#	    using CMake.
#
#
#	2.) This setup does not actually compile a library.  Instead, the appropriate
#	    files are added to PARENT_SCOPE variables so that an upper-level
#	    CMakeLists.txt file may properly build.
#





# ==============================================================================
#
#   Set basic options
#
# ==============================================================================

# Require CMake version 2.8+
cmake_minimum_required(VERSION 2.8)

# Set name ofproject
set(GENCLASSES_NAME "genCppClasses" PARENT_SCOPE)


# Control build options and feature support.  Initial values are inherited from
# parent directory.
SET(GENCLASSES_CXX_11 ${USE_CXX_11} CACHE BOOL "Enable C++-11 Support within General-Use C++ classes" FORCE)
SET(GENCLASSES_OPENMP ${USE_OPENMP} CACHE BOOL "Enable OpenMP Support within General-Use C++ classes" FORCE)
SET(GENCLASSES_USE_QT ${USE_QT} CACHE BOOL "Enable Qt Functionality within General-Use C++ classes" FORCE)
SET(GENCLASSES_USE_VTK ${USE_VTK} CACHE BOOL "Enable VTK Libraries within General-Use C++ classes" FORCE)

# Mark these variables as advanced since they likely won't need to be adjusted
MARK_AS_ADVANCED(GENCLASSES_CXX_11 GENCLASSES_OPENMP GENCLASSES_USE_QT GENCLASSES_USE_VTK)

if(GENCLASSES_USE_VTK)
    MARK_AS_ADVANCED(CLEAR VTK_DIR)
else()
    MARK_AS_ADVANCED(VTK_DIR)
endif()    







# ==============================================================================
#
#   Set compiler flags
#
# ==============================================================================


# GCC
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if(GENCLASSES_CXX_11)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")    # Enable C++11 features
        add_definitions(-DCXX11)
    endif()

    if(GENCLASSES_OPENMP)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")      # Enable OpenMP
    endif()
endif()


# MSVC
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    if(GENCLASSES_CXX_11)
        add_definitions(-DCXX11)
    endif()

    if(GENCLASSES_OPENMP)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")      # Enable OpenMP
    endif()
endif()






# ==============================================================================
#
#   Define header/source files
#
# ==============================================================================


# Define source files to be compiled into object files.  
set(GENCLASSES_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/ProgressInfo.cxx
                      ${CMAKE_CURRENT_SOURCE_DIR}/QtIntermediaryBase.cxx
                      ${CMAKE_CURRENT_SOURCE_DIR}/QtIntermediary.cxx
                      ${CMAKE_CURRENT_SOURCE_DIR}/StringManip.cxx
    )



# Define header files
include_directories(${CMAKE_CURRENT_SOURCE_DIR})   # Include source directory during compilation
file(GLOB_RECURSE GENCLASSES_HEADERS "*.h")        # Add all *.h files to "GEN_CLASSES_HEADERS"
file(GLOB_RECURSE GENCLASSES_CPP "*.cpp")          # Add all *.cpp files to "GEN_CLASSES_CPP"


# Add header and source files to project-wide variables.
set(HEADERS ${HEADERS} ${GENCLASSES_HEADERS} PARENT_SCOPE)
set(SOURCES ${SOURCES} ${GENCLASSES_SOURCE} PARENT_SCOPE)

# NOTE: Files added to 'SOURCES' are 'compileable' sources, not templated sources which are
#	included in their header files.







# ==============================================================================
#
#   Set variables based on configuration options
#
# ==============================================================================

# Use Qt
if(GENCLASSES_USE_QT)
    # Require Qt4
    FIND_PACKAGE(Qt4 REQUIRED)
    INCLUDE(${QT_USE_FILE})
    ADD_DEFINITIONS(${QT_DEFINITIONS})

    SET(QT_HEADERS ${QT_HEADERS}
            QSafeApplication.h
            QtIntermediary.h
            QtIntermediaryBase.h
            PARENT_SCOPE
       )
    # These headers are not actually processed yet.  This CMakeLists.txt file only
    # adds them to 'QT_HEADERS' so that they can be processed at a later time

endif()

# Use VTK
if(GENCLASSES_USE_VTK)
    FIND_PACKAGE(VTK REQUIRED)
    INCLUDE(${VTK_USE_FILE})
    
    add_definitions(-DUSE_VTK)

endif()
