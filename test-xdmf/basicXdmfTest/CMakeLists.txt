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

add_library(basicXdmfTest STATIC basicXdmfTest.cpp)
target_link_libraries(basicXdmfTest ${XDMF_LIBS})

set(TESTVAR "test")
message(STATUS "XDMF Libs: "${XDMF_LIBS})
