# ==============================================================================
#
#   Set basic options
#
# ==============================================================================

# Require CMake version 2.8+
cmake_minimum_required(VERSION 2.8)

# Set name ofproject
set(GEN_CLASSES_NAME "genCppClasses")

# Set default build type
set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)


# Control build options and feature support
SET(GEN_CLASSES_CXX_11 ON CACHE BOOL "Enable C++-11 Support")
SET(GEN_CLASSES_OPENMP ON CACHE BOOL "Enable OpenMP Support")
SET(GEN_CLASSES_USE_QT ON CACHE BOOL "Enable Qt Functionality")
SET(GEN_CLASSES_USE_VTK ON CACHE BOOL "Enable VTK Libraries")








# ==============================================================================
#
#   Set compiler flags
#
# ==============================================================================


# GCC
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if(GEN_CLASSES_CXX_11)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")    # Enable C++11 features
        add_definitions(-DCXX11)
    endif()

    if(GEN_CLASSES_OPENMP)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")      # Enable OpenMP
    endif()
endif()


# MSVC
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    if(GEN_CLASSES_CXX_11)
        add_definitions(-DCXX11)
    endif()

    if(GEN_CLASSES_OPENMP)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")      # Enable OpenMP
    endif()
endif()






# ==============================================================================
#
#   Define header/source files
#
# ==============================================================================


# Define source files to be compiled into object files.  Most of the *.cpp files are not
# listed here because they are templated and thus directly included in their headers.
set(GEN_CLASSES_LIB_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/ProgressInfo.cpp
                           ${CMAKE_CURRENT_SOURCE_DIR}/QtIntermediaryBase.cpp
                           ${CMAKE_CURRENT_SOURCE_DIR}/QtIntermediary.cpp
                           ${CMAKE_CURRENT_SOURCE_DIR}/StringManip.cpp
    )



# Define header files
include_directories(${CMAKE_CURRENT_SOURCE_DIR})    # Include source directory during compilation
file(GLOB_RECURSE GEN_CLASSES_HEADERS "*.h")        # Add all *.h files to "GEN_CLASSES_HEADERS"
file(GLOB_RECURSE GEN_CLASSES_CPP "*.cpp")          # Add all *.cpp files to "GEN_CLASSES_CPP"





# ==============================================================================
#
#   Specify library to be compiled & related source files
#
# ==============================================================================

add_library(${GEN_CLASSES_NAME} SHARED
                    ${GEN_CLASSES_LIB_SOURCE} ${GEN_CLASSES_HEADERS}
                    ${QT_HEADERS}
           )





# ==============================================================================
#
#   Set variables based on configuration options
#
# ==============================================================================

# Use Qt
if(GEN_CLASSES_USE_QT)
    # Require Qt4
    FIND_PACKAGE(Qt4 REQUIRED)
    INCLUDE(${QT_USE_FILE})
    ADD_DEFINITIONS(${QT_DEFINITIONS})

    SET(QT_HEADERS
            QSafeApplication.h
            QtIntermediary.h
            QtIntermediaryBase.h
       )
    QT4_WRAP_CPP(QT_HEADERS_MOC ${QT_HEADERS})

    if(QT_LIBRARIES)
      target_link_libraries(${GEN_CLASSES_NAME} ${QT_LIBRARIES})
    endif()
endif()

# Use VTK
if(GEN_CLASSES_USE_VTK)
    FIND_PACKAGE(VTK REQUIRED)
    INCLUDE(${VTK_USE_FILE})

    if(VTK_LIBRARIES)
      target_link_libraries(${GEN_CLASSES_NAME} ${VTK_LIBRARIES})
    else()
      target_link_libraries(${GEN_CLASSES_NAME} vtkHybrid)
    endif()
endif()





# ==============================================================================
#
#   Set install rules
#
# ==============================================================================

# Move header files to "include" within install directory.  The *.cpp files
# are also included here since most of them are included in the header files.
# Strictly speaking, only those included in headers need to be copied.  The
# *.cpp files which are compiled do not need to be moved.
INSTALL(FILES ${GEN_CLASSES_HEADERS} ${GEN_CLASSES_CPP} DESTINATION include)

# Move runtime objects to "bin", and library/archive objects to "lib" within
# the install directory.
INSTALL(TARGETS ${GEN_CLASSES_NAME}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
