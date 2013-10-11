# ==============================================================================
#
#   NOTES
#
# ==============================================================================
#
#   - The purpose of this file is to configure CMake to properly handle the
#     code within this directory.  It is expected that this code is being
#     used as part of a larger project and is not the project itself.
#
#   - To incorporate the code in this directory into a project, there are two
#     options:
#       1.) Copy this source directory to be a subdirectory of the project
#           source and use the 1-argument version of 'ADD_SUBDIRECTORY()'.
#           This approach will make integration with version control either
#           difficult or duplicative due to the copying of this code.
#       2.) Add this directory to the project using the 2-argument version
#           of 'ADD_SUBDIRECTORY()'.  This approach is preferred since it
#           requires no copying of source code and does not complicate
#           the management of the source code repository.
#
#   - Several configuration options are inherited from the project itself,
#     such as the use of C++11, VTK support, OpenMP support, etc.  Enabling
#     support for such features within this directory's code requires *first*
#     enabling the desired features at the project-level.  For example, if
#     OpenMP is disabled project-wide, it will not be supported by the code
#     in this directory and that cannot be overridden.  However, if OpenMP
#     is enabled in the project, it does not necessarily have to be enabled
#     here.
#
#   - Several variables are exported to the parent directory using 'PARENT_SCOPE'.
#     These variables are then accessible within the parent directory:
#
#           HEADERS             Standard C/C++ header files
#
#           HEADERS_QT          Qt-derived header files which require MOC
#
#           RESOURCES_QT        Qt-based resource files which require MOC
#
#           SOURCES             Standard C/C++ source files which require compilation
#
#           SOURCES_CUDA        CUDA source files (e.g., *.cu) which are to be
#                               included into a single project-wide *.cu file
#                               for NVCC to compile
#
#           SOURCES_QT          Qt-based C++ source files.  These do not require
#                               MOC and are ready for normal compilation along
#                               with 'SOURCES'
#
#           SOURCES_TEMPLATE    Standard C/C++ source files which are *NOT* to be
#                               compiled due to their containing of templated
#                               functions.  These files are expected to be
#                               #include'd at the end of their respective headers
#
#           UI_QT               Qt-based interface-definition files (e.g., *.ui)
#                               which require MOC
#
#   - To configure this CMakeLists.txt template for a source directory, do the following:
#       1.) Rename this file to "CMakeLists.txt"
#       2.) Set variables in the first two sections (basic options & file lists)
#           - If no files are present for a given category, simply comment-out
#             the lines used for the SET() command.  Comment symbol is hash: "#".
#
#
#




# ==============================================================================
#
#   Set basic options
#
# ==============================================================================

# Set name of code in this directory.  NO SPACES!
SET(CODE_NAME "GenClasses")

# Do a Find/Replace in this entire file:
#   Find: GenClasses
#   Replace: CodeName (Replace with the code name above)
# Also update names in the generated header input file.


# Set human-readable name of code in this directory.  Spaces are OK.
SET(CODE_NAME_HUMAN "General C++ Classes")

# Do a Find/Replace in this entire file:
#   Find: General C++ Classes
#   Replace: CodeName (Replace with the human-readable code name above)


# Set name of developer for code in this directory.
SET(CODE_DEVELOPER "RobertGrandin")


# Set name of input/output for generated header file
SET(GENERATED_HEADER_INPUT "GenClassesConfig.h.in")
SET(GENERATED_HEADER_OUTPUT "GenClassesConfig.h")




# ==============================================================================
#
#   Define header/source files for PARENT_SCOPE variables
#
# ==============================================================================


# Standard headers -- These require no extra processing.  Includes headers for
#                     templated classes.
include_directories(${CMAKE_CURRENT_SOURCE_DIR}) # Include source directory during compilation.
set(HEADERS ${HEADERS}
            PARENT_SCOPE)


# Qt headers -- These require MOC pre-processing before compilation.
set(HEADERS_QT ${HEADERS_QT}
               "${CMAKE_CURRENT_LIST_DIR}/QSafeApplication.h"
               "${CMAKE_CURRENT_LIST_DIR}/QtIntermediary.h"
               "${CMAKE_CURRENT_LIST_DIR}/QtIntermediaryBase.h"
               PARENT_SCOPE)


# Qt resources --  These require MOC pre-processing before compilation.
set(RESOURCES_QT ${RESOURCES_QT}
                 PARENT_SCOPE)


# Standard sources --  These require no extra processing.
set(SOURCES ${SOURCES}
            "${CMAKE_CURRENT_LIST_DIR}/ProgressInfo.cxx"
            "${CMAKE_CURRENT_LIST_DIR}/StringManip.cxx"
            PARENT_SCOPE)


# CUDA sources -- These require the NVCC compiler.
set(SOURCES_CUDA ${SOURCES_CUDA}
                 PARENT_SCOPE)


# Qt sources --  These require no extra processing.
set(SOURCES_QT ${SOURCES_QT}
               "${CMAKE_CURRENT_LIST_DIR}/QtIntermediary.cxx"
               "${CMAKE_CURRENT_LIST_DIR}/QtIntermediaryBase.cxx"
               PARENT_SCOPE)


# Templated sources -- These are NOT to be compiled.  They should be #include'd at the
#                      end of their respective header files.  This variable exists to
#                      provide a convenient list/reference.
set(SOURCES_TEMPLATE ${SOURCES_TEMPLATE}
                     PARENT_SCOPE)


# Qt user-interface definition files --  These require MOC pre-processing before compilation.
set(UI_QT ${UI_QT}
          PARENT_SCOPE)










# ==============================================================================
#
#   NO EDITS REQUIRED BELOW HERE
#
# ==============================================================================

# Require CMake version 2.8+
cmake_minimum_required(VERSION 2.8)






# ==============================================================================
#
#   Set local enable/disable options for feature-support & check against
#   global option
#
# ==============================================================================

# Set values.  Default to global settings.
SET(${GenClasses_ENABLE_CXX11} ${USE_CXX_11} CACHE BOOL "Enable C++-11 support within General C++ Classes")
SET(${GenClasses_ENABLE_OPENMP} ${USE_OPENMP} CACHE BOOL "OpenMP support within General C++ Classes")
SET(${GenClasses_ENABLE_VTK} ${USE_VTK} CACHE BOOL "Enable VTK support within General C++ Classes")
SET(${GenClasses_ENABLE_CUDA} ${USE_CUDA} CACHE BOOL "Enable CUDA support within General C++ Classes")
SET(${GenClasses_ENABLE_QT} ${USE_QT} CACHE BOOL "Enable QT support within General C++ Classes")


# Mark these as advanced options.
MARK_AS_ADVANCED(FORCE
                 ${GenClasses_ENABLE_CXX11}
                 ${GenClasses_ENABLE_OPENMP}
                 ${GenClasses_ENABLE_VTK}
                 ${GenClasses_ENABLE_CUDA}
                 ${GenClasses_ENABLE_QT}
                )


# Check that local option is compatible with global option.  If a feature is globally
# disabled it cannot be enabled here.  If a feature is globally enabled, it does not
# have to be enabled here.
IF(NOT ${USE_CXX_11})
    SET(${GenClasses_ENABLE_CXX11} FALSE)
ENDIF()

IF(NOT ${USE_OPENMP})
    SET(${GenClasses_ENABLE_OPENMP} FALSE)
ENDIF()

IF(NOT ${USE_VTK})
    SET(${GenClasses_ENABLE_VTK} FALSE)
ENDIF()

IF(NOT ${USE_CUDA})
    SET(${GenClasses_ENABLE_CUDA} FALSE)
ENDIF()

IF(NOT ${USE_QT})
    SET(${GenClasses_ENABLE_QT} FALSE)
ENDIF()






# ==============================================================================
#
#   Setup CMake-configured header file
#
# ==============================================================================

CONFIGURE_FILE(
    "${CMAKE_CURRENT_LIST_DIR}/${GENERATED_HEADER_INPUT}"
    "${CMAKE_CURRENT_LIST_DIR}/${GENERATED_HEADER_OUTPUT}"
    )
