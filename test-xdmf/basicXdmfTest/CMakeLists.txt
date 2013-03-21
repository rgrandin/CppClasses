set(basicXdmfTestName "basicXdmfTest")



set(LIB_SOURCE ${GEN_CXX_CLASSES_DIR}/QtIntermediaryBase.cpp
               ${GEN_CXX_CLASSES_DIR}/StringManip.cpp
    )


# Add project dependencies -- XDMF library
find_package(XDMF REQUIRED)
if(${XDMF_FOUND})
    #include_directories(${XDMF_INCLUDE_PATH})
    #target_link_libraries(${PROJECT_NAME} ${XDMF_LIBRARIES})
    set(XDMF_INCLUDES ${XDMF_INCLUDE_PATH})
    set(XDMF_LIBS ${XDMF_LIBRARIES})
else()
    set(XDMF_INCLUDES "/usr/local/include")
    set(XDMF_LIBS "-L/usr/local/xdmf/lib/ -lmetis -lvtkexoIIc -lvtkhdf5 -lvtklibxml2 -lXdmf -lXdmfUtils")
endif()

# Manually set includes and lib variables.  This is because the 'find_package' command doesn't seem to
# work for Xdmf right now.
set(XDMF_INCLUDES "/usr/local/include")
set(XDMF_LIBS "-L/usr/local/xdmf/lib/ -lmetis -lvtkexoIIc -lvtkhdf5 -lvtklibxml2 -lXdmf -lXdmfUtils")

include_directories(${XDMF_INCLUDES})

message(STATUS "PROJECT_SOURCE_DIR: "${PROJECT_SOURCE_DIR})
message(STATUS "GEN_CXX_CLASSES_DIR: "${GEN_CXX_CLASSES_DIR})
message(STATUS "LIB_SOURCE: "${LIB_SOURCE})

add_library(basicXdmfTest SHARED basicXdmfTest.cpp ${LIB_SOURCE})
target_link_libraries(basicXdmfTest ${XDMF_LIBS})

#set(TESTVAR "test")
#message(STATUS "XDMF Libs: "${XDMF_LIBS})



# Add project dependencies -- HDF5 library
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_PATH})

# Enable C++11 support for GCC compiler
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()




