set(basicXdmfTestName "basicXdmfTest")          # Define library  name



# Add project dependencies -- XDMF library
find_package(XDMF REQUIRED)
if(${XDMF_FOUND})
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

# Add project dependencies -- HDF5 library
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_PATH})

include_directories(${XDMF_INCLUDES})
include_directories(${GEN_CPP_CLASSES_DIR})


#message(STATUS "PROJECT_SOURCE_DIR: "${PROJECT_SOURCE_DIR})
#message(STATUS "GEN_CXX_CLASSES_DIR: "${GEN_CPP_CLASSES_DIR})


add_library(${basicXdmfTestName} SHARED basicXdmfTest.cpp)
target_link_libraries(${basicXdmfTestName} ${XDMF_LIBS})

include_directories(${GEN_CPP_CLASSES_DIR})
add_subdirectory(${GEN_CPP_CLASSES_DIR} "${PROJECT_BINARY_DIR}/gen-cpp-classes")
target_link_libraries(${basicXdmfTestName} "genCppClasses")

