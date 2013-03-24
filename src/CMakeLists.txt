set(GEN_CLASSES_NAME "genCppClasses")                       # Define name of library


# Set GCC compiler flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")    # Enable C++11 features
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")      # Enable OpenMP
endif()



set(GEN_CPP_CLASSES_DIR ${CMAKE_CURRENT_SOURCE_DIR})        # Define directory containing source files to be the
#message(STATUS "General C++ Dir: "${GEN_CPP_CLASSES_DIR})  # same directory as this CMakeLists.txt file


# Define source files to be compiled into object files.  Most of the *.cpp files are not
# listed here because they are templated and thus directly included in their headers.
set(GEN_CLASSES_LIB_SOURCE ${GEN_CPP_CLASSES_DIR}/QtIntermediaryBase.cpp
                           ${GEN_CPP_CLASSES_DIR}/QtIntermediary.cpp
                           ${GEN_CPP_CLASSES_DIR}/StringManip.cpp
    )




# Add library for general C++ classes
add_library(${GEN_CLASSES_NAME} SHARED ${GEN_CLASSES_LIB_SOURCE})
