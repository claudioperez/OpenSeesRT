#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# Claudio Perez
#==============================================================================

cmake_minimum_required(VERSION 3.16)
project(OpenSees C CXX Fortran)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(OPS_EXTERN_SOURCE_DIR "${PROJECT_SOURCE_DIR}/OTHER/")
set(OPS_EXTERNALS_DIR "${PROJECT_SOURCE_DIR}/OTHER/")

set(OPS_BUNDLED_DIR "${PROJECT_SOURCE_DIR}/OTHER/")
set(OPS_SRC_DIR "${PROJECT_SOURCE_DIR}/SRC/")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/BuildTools/cmake")
include(OpenSeesFunctions)
enable_testing()

add_definitions(-DOPS_USE_RUNTIME -D_THREADS)

# add Conan stuff for external dependencies
if(EXISTS "${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
   message(">>> Using Conan")
   set(USE_CONAN TRUE)
   include("${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
   conan_basic_setup(TARGETS)
   set(TCL_LIBRARIES "CONAN_PKG::tcl")
   set(TCL_LIBRARY "CONAN_PKG::tcl")
endif()

# pseudo target modeling all exteral packages
add_library(OPS_External_packages INTERFACE)

# include user config
include(${PROJECT_SOURCE_DIR}/Conf.cmake)

# define user-selectable options for end target
set_property(CACHE OPS_FINAL_TARGET PROPERTY STRINGS 
    G3 OpenSeesTcl OpenSeesMP OpenSeesSP)

#==============================================================================
#                            OS Configuration
#==============================================================================
add_library(OPS_OS_Specific_libs INTERFACE)

if (UNIX)
  if (DEFINED Dependencies)
    include("OpenSeesDependencies${Dependencies}")
  else ()
    include(OpenSeesDependenciesUnix)
  endif ()

  if (APPLE)
   message(STATUS ">>> MacOS")
   add_compile_definitions(_LINUX _UNIX _TCL85)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmax-errors=1")
  else()
     message(STATUS ">>> LINUX")
     add_compile_definitions(_LINUX _UNIX USE_TCL_STUBS _TCL85)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmax-errors=1 -ffloat-store")
  endif()

  add_compile_options(
    $<$<CONFIG:DEBUG>:-ggdb>
    $<$<CONFIG:DEBUG>:-g3>
    $<$<CONFIG:DEBUG>:-Og>
    $<$<CONFIG:RELEASE>:-O3>
  )
endif()

if(WIN32) # NOTE: this will execute for both 32-bit and 64-bit builds.
    include(OpenSeesDependenciesWin)
    # target_link_libraries(OPS_OS_Specific_libs INTERFACE wsock32 ws2_32)
    add_compile_definitions(_WIN32 _TCL85)
    message(STATUS ">>> WIN32")
endif()

message("OPS >>> BLAS:    ${BLAS_LIBRARIES}")
message("OPS >>> CBLAS:   ${CBLAS_LIBRARY}\n")
message("OPS >>> LAPACK:  ${LAPACK_LIBRARIES}")
message("OPS >>> SUPERLU: ${SUPERLU_LIBRARIES}")
message("OPS >>> ARPACK:  ${ARPACK_LIBRARIES}")
message("OPS >>> UMFPACK: ${UMFPACK_LIBRARIES}")
message("OPS >>> CSPARSE: ${CSPARSE_LIBRARIES}")
message("OPS >>> TCL:     ${TCL_LIBRARIES}")
message("                 ${TCL_LIBRARY}  ")
message("                $ENV{TCL_LIBRARY}")

message("OPS >>> AMD:     ${AMD_LIBRARIES}")


#==============================================================================
#                            General Setup
#
#==============================================================================
#----------------------------------------------------------------
# Compilers
#---------------------------------------------------------------- 

# Fortran
#--------------------------------------
enable_language(Fortran)

# C++
#--------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> cqls <TARGET> <LINK_FLAGS> <OBJECTS>")
add_compile_options(
  # "$<$<CXX_COMPILER_ID:MSVC>:/W4>"
  "$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall>"
  "$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wextra>"
  "$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-unused-parameter>"
  "$<$<COMPILE_LANG_AND_ID:CXX,AppleClang,Clang>:-Wno-reorder>"
  "$<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wno-reorder>"
  # "$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wpedantic>"
  # $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wsign-conversion>
  # $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wconversion>
)

#----------------------------------------------------------------
# Global Includes
#---------------------------------------------------------------- 
# TODO BeamIntegration and MatrixUtil need to be removed from element/forceBeamColumn

include_directories(
  ${OPS_SRC_DIR}
  ${OPS_SRC_DIR}/actor/actor
  ${OPS_SRC_DIR}/actor/objectBroker
  ${OPS_SRC_DIR}/actor/channel
  ${OPS_SRC_DIR}/actor/message

  ${OPS_SRC_DIR}/analysis/algorithm
  ${OPS_SRC_DIR}/analysis/algorithm/domainDecompAlgo
  ${OPS_SRC_DIR}/analysis/algorithm/eigenAlgo
  ${OPS_SRC_DIR}/analysis/algorithm/equiSolnAlgo
  ${OPS_SRC_DIR}/analysis/algorithm/equiSolnAlgo/accelerator
  ${OPS_SRC_DIR}/analysis/analysis
  ${OPS_SRC_DIR}/analysis/dof_grp
  ${OPS_SRC_DIR}/analysis/fe_ele
  ${OPS_SRC_DIR}/analysis/handler
  ${OPS_SRC_DIR}/analysis/integrator
  ${OPS_SRC_DIR}/analysis/model
  ${OPS_SRC_DIR}/analysis/model/simple
  ${OPS_SRC_DIR}/analysis/numberer

  ${OPS_SRC_DIR}/interpreter
  ${OPS_SRC_DIR}/api
  ${OPS_SRC_DIR}/convergenceTest
  ${OPS_SRC_DIR}/coordTransformation
  ${OPS_SRC_DIR}/database
  ${OPS_SRC_DIR}/domain/component
  ${OPS_SRC_DIR}/domain/constraints
  ${OPS_SRC_DIR}/domain/domain
  ${OPS_SRC_DIR}/domain/domain/partitioned
  ${OPS_SRC_DIR}/domain/domain/single
  ${OPS_SRC_DIR}/domain/groundMotion
  ${OPS_SRC_DIR}/domain/load
  ${OPS_SRC_DIR}/domain/loadBalancer
  ${OPS_SRC_DIR}/domain/node
  ${OPS_SRC_DIR}/domain/partitioner
  ${OPS_SRC_DIR}/domain/pattern
  ${OPS_SRC_DIR}/domain/pattern/drm
  ${OPS_SRC_DIR}/domain/region
  ${OPS_SRC_DIR}/domain/subdomain
  ${OPS_SRC_DIR}/element
  ${OPS_SRC_DIR}/graph/graph
  ${OPS_SRC_DIR}/graph/numberer
  ${OPS_SRC_DIR}/graph/partitioner
  ${OPS_SRC_DIR}/handler
  ${OPS_SRC_DIR}/material
  ${OPS_SRC_DIR}/material/nD
  ${OPS_SRC_DIR}/material/Solid
  ${OPS_SRC_DIR}/material/section
  ${OPS_SRC_DIR}/material/section/fiber
  ${OPS_SRC_DIR}/material/section/yieldSurface
  ${OPS_SRC_DIR}/material/uniaxial
  ${OPS_SRC_DIR}/matrix
  ${OPS_SRC_DIR}/modelbuilder
  ${OPS_SRC_DIR}/modelbuilder/tcl
  ${OPS_SRC_DIR}/recorder
  ${OPS_SRC_DIR}/recorder/response
  ${OPS_SRC_DIR}/system_of_eqn
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE
  ${OPS_SRC_DIR}/system_of_eqn/eigenSOE
  ${OPS_SRC_DIR}/tagged
  ${OPS_SRC_DIR}/tagged/storage
  ${OPS_SRC_DIR}/utility
  ${OPS_BUNDLED_DIR}/SuperLU_5.1.1/SRC/
)


include_directories(${TCL_INCLUDE_PATH})


#==============================================================================
#                            Define Targets
#
#==============================================================================
add_library(OPS_Numerics           INTERFACE)
target_link_libraries(OPS_Numerics INTERFACE
  ${ARPACK_LIBRARIES} 
# ${SUPERLU_LIBRARIES}
  ${AMD_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES} 
)

# Core OpenSees
add_library(OPS_Matrix             OBJECT)
add_library(OPS_Actor              OBJECT)
add_library(OPS_ObjectBroker       OBJECT)
add_library(G3_ObjectBroker        OBJECT)
add_library(OPS_Handler            OBJECT)
add_library(OPS_Recorder           OBJECT)
add_library(OPS_Tagged             OBJECT)
add_library(OPS_Utilities          OBJECT)
add_library(OPS_ModelBuilder       OBJECT)
add_library(OPS_Domain             OBJECT)
add_library(OPS_SysOfEqn           OBJECT)

add_library(OPS_Analysis           OBJECT)
add_library(OPS_ConvergenceTest    OBJECT)
add_library(OPS_Thermal            OBJECT)
add_library(OPS_Element            OBJECT)
add_library(OPS_ElementFortran     OBJECT)
add_library(OPS_Material           OBJECT)
add_library(OPS_MaterialFortran    OBJECT)
add_library(OPS_Damage             OBJECT)
add_library(OPS_Database           OBJECT)

# Optional Extensions
add_library(OPS_ASDEA              OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Paraview           OBJECT EXCLUDE_FROM_ALL)


opensees_library(OPS_SysOfEqn_UMF
  LINK ${UMFPACK_LIBRARIES} ${AMD_LIBRARIES})

# Packaged libraries
add_library(G3) 

set_property(TARGET OPS_MaterialFortran   PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET OPS_ElementFortran    PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET G3                    PROPERTY POSITION_INDEPENDENT_CODE 1)


# Executables and final targets
#------------------------------------------------------------------------------
# Adding EXCLUDE_FROM_ALL prevents these from being made by default.

add_executable(OpenSeesSP    EXCLUDE_FROM_ALL )
add_executable(OpenSeesMP    EXCLUDE_FROM_ALL )


# Add sources to targets
add_subdirectory(${OPS_SRC_DIR})

#==============================================================================
#                            Configure targets
#
#==============================================================================
target_link_libraries(OPS_Numerics     INTERFACE ${OPS_Numlib_List})
message(">>> Numeric libs: ${OPS_Numlib_List}")
target_link_libraries(OPS_Element      PRIVATE ${OPS_Element_List})
target_link_libraries(OPS_ObjectBroker PRIVATE ${OPS_Element_List})

target_link_libraries(OPS_Material PRIVATE ${TCL_LIBRARY})


# G3
#----------------------------
target_link_libraries(G3 
    OPS_Database
    OPS_Matrix
    OPS_Actor
    OPS_Analysis 
    OPS_Domain
    OPS_Thermal
    OPS_ConvergenceTest
    OPS_Element
    OPS_ElementFortran
    OPS_Material
    OPS_MaterialFortran
    OPS_Material_YieldSurface
    OPS_Section_Repres
    OPS_Section_YieldSurface
    OPS_Recorder
    OPS_Handler
    OPS_SysOfEqn
    ${OPS_SysOfEqn_List}
    OPS_Tagged
    OPS_Utilities
    graph
    coordTransformation
    G3_ObjectBroker
    OPS_Numerics 
)


#==============================================================================
#                            Apply Options
#
#==============================================================================

# Extensions
#----------------------------
message("OPS >>> Configuring OpenSees extensions")
foreach(extension IN LISTS OPS_SysOfEqn_List OPS_Element_List OPS_Extension_List)
    string(TOUPPER "${extension}" ext_flag) 
    string(REGEX REPLACE "^OPS_" "OPSDEF_" ext_flag "${ext_flag}")
    add_compile_definitions(${ext_flag})
endforeach()

foreach(extension IN LISTS OPS_Exclude_List)
    string(TOUPPER "${extension}" ext_flag) 
    string(REGEX REPLACE "^OPS_" "OPS_EXCLUDE_" ext_flag "${ext_flag}")
    message("    Adding macro definition '${ext_flag}'")
    add_compile_definitions(${ext_flag})
endforeach()


add_compile_definitions(_NOGRAPHICS)


# Reliability
#----------------------------
if ("OPS_Reliability" IN_LIST OPS_Extension_List)
    add_compile_definitions(_RELIABILITY)
    target_link_libraries(${OPS_FINAL_TARGET} OPS_Reliability)
endif()


# HDF5
#----------------------------
if(OPS_Use_HDF5)
   find_package(HDF5)
    if(HDF5_FOUND)
        include_directories(${HDF5_INCLUDE_DIR})
        set(_hdf5_libs hdf5 hdf5_cpp)
    add_compile_definitions(-D_H5DRM)
    else()
        message(STATUS "OPS >>> Could not find HDF5")
    endif()
endif()


# Include test suite
#add_subdirectory(EXAMPLES/)


get_target_property(OPS_Material_COMPILE_OPTIONS OPS_Material COMPILE_OPTIONS)
  string(REPLACE "-Wall" "" OPS_Material_COMPILE_OPTIONS "${OPS_Material_COMPILE_OPTIONS}")
  string(REPLACE "-Wextra" "" OPS_Material_COMPILE_OPTIONS "${OPS_Material_COMPILE_OPTIONS}")
  set_target_properties(OPS_Material PROPERTIES COMPILE_OPTIONS "${OPS_Material_COMPILE_OPTIONS}")

get_target_property(OPS_Element_COMPILE_OPTIONS OPS_Element COMPILE_OPTIONS)
  string(REPLACE "-Wall" "" OPS_Element_COMPILE_OPTIONS "${OPS_Element_COMPILE_OPTIONS}")
# string(REPLACE "-Wextra" "" OPS_Element_COMPILE_OPTIONS "${OPS_Element_COMPILE_OPTIONS}")
  set_target_properties(OPS_Element PROPERTIES COMPILE_OPTIONS "${OPS_Element_COMPILE_OPTIONS}")

