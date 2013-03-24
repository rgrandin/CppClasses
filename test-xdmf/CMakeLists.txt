# Require CMake version 2.8+
cmake_minimum_required(VERSION 2.8)


# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()


# Set project info
project(xdmf-demo)                  # Set name of project
set(XDMF_DEMO_VERSION_MAJOR 1)      # Set major version number
set(XDMF_DEMO_VERSION_MINOR 0)      # Set minor version number


aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})



# Set file to contain version information
configure_file(
  "${PROJECT_SOURCE_DIR}/XdmfDemoConfig.h.in"
  "${PROJECT_BINARY_DIR}/XdmfDemoConfig.h"
  )
include_directories("${PROJECT_BINARY_DIR}")    # Include req'd to find 'XdmfDemoConfig.h'


# Set GCC compiler flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")    # Enable C++11 features
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")      # Enable OpenMP
endif()


# Set build type
#set(CMAKE_BUILD_TYPE "RelWithDebInfo")
set(CMAKE_BUILD_TYPE "Debug")


# Add project dependencies -- MPI
find_package(MPI REQUIRED)
if(MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARY} ${MPI_CXX_LIBRARIES})
endif()


# Add General-Use C++ Classes
set(GEN_CXX_CLASSES_DIR "${PROJECT_SOURCE_DIR}/../src")     # Manually define directory
include_directories(${GEN_CXX_CLASSES_DIR})





# Add project dependencies -- My basic test of XDMF library
include_directories("${PROJECT_SOURCE_DIR}/basicXdmfTest")
add_subdirectory("basicXdmfTest")
target_link_libraries(${PROJECT_NAME} basicXdmfTest)

