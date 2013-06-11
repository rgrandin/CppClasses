cmake_minimum_required(VERSION 2.8)


set(GEN_CLASSES_NAME "genCppClasses")                       # Define name of library


set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)


SET(GEN_CLASSES_CXX_11 ON CACHE BOOL "Enable C++-11 Support")
SET(GEN_CLASSES_OPENMP ON CACHE BOOL "Enable OpenMP Support")


# Set GCC compiler flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if(GEN_CLASSES_CXX_11)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")    # Enable C++11 features
        add_definitions(-DCXX11)
    endif()

    if(GEN_CLASSES_OPENMP)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")      # Enable OpenMP
    endif()
endif()


# Set MSVC compiler flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    if(GEN_CLASSES_CXX_11)
        add_definitions(-DCXX11)
    endif()

    if(GEN_CLASSES_OPENMP)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")      # Enable OpenMP
    endif()
endif()


set(GEN_CPP_CLASSES_DIR ${CMAKE_CURRENT_SOURCE_DIR})        # Define directory containing source files to be the



# Define source files to be compiled into object files.  Most of the *.cpp files are not
# listed here because they are templated and thus directly included in their headers.
set(GEN_CLASSES_LIB_SOURCE ${GEN_CPP_CLASSES_DIR}/ProgressInfo.cpp
                           ${GEN_CPP_CLASSES_DIR}/QtIntermediaryBase.cpp
                           ${GEN_CPP_CLASSES_DIR}/QtIntermediary.cpp
                           ${GEN_CPP_CLASSES_DIR}/StringManip.cpp
    )


# Define header files
include_directories(${GEN_CPP_CLASSES_DIR})
file(GLOB_RECURSE GEN_CLASSES_HEADERS "*.h")



# Add library for general C++ classes
add_library(${GEN_CLASSES_NAME} SHARED ${GEN_CLASSES_LIB_SOURCE} ${GEN_CLASSES_HEADERS})
