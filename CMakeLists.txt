#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# Claudio Perez
#==============================================================================

cmake_minimum_required(VERSION 3.16)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(OpenSees C CXX Fortran)

set(OPS_BUNDLED_DIR   "${CMAKE_CURRENT_LIST_DIR}/OTHER/")
set(OPS_SRC_DIR       "${CMAKE_CURRENT_LIST_DIR}/SRC/")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}"
                      "${CMAKE_CURRENT_LIST_DIR}/BuildTools/cmake")

include(OpenSeesFunctions)

enable_testing()

add_compile_definitions(OPS_USE_RUNTIME USE_TCL_STUBS _TCL85 _NOGRAPHICS) # _THREADS

# add Conan stuff for external dependencies
if(EXISTS "${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
   message(":: Using Conan")
   set(USE_CONAN TRUE)
   include("${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
   conan_basic_setup(TARGETS)
   set(TCL_INCLUDE_DIRS "${CONAN_INCLUDE_DIRS_TCL}")
   set(TCL_INCLUDE_PATH "${CONAN_INCLUDE_DIRS_TCL}")
   set(TCL_LIBRARIES    "${CONAN_PKG_LIBS_TCL}")

elseif(WIN32)
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})
  list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR})
  
  if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
    message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
    file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/0.18.1/conan.cmake"
                  "${CMAKE_BINARY_DIR}/conan.cmake"
                  TLS_VERIFY ON)
  endif()
  
  include(${CMAKE_BINARY_DIR}/conan.cmake)
  
  conan_cmake_configure(REQUIRES tcl/8.6.11
                        GENERATORS cmake ) # cmake_find_package)

  conan_cmake_autodetect(settings)
  
  conan_cmake_install(PATH_OR_REFERENCE ${CMAKE_SOURCE_DIR}
                      BUILD missing
                      REMOTE conancenter
                      SETTINGS ${settings})

  include(${CMAKE_BINARY_DIR}/src/libg3/conanbuildinfo.cmake)
  conan_basic_setup() # TARGETS)


  set(TCL_INCLUDE_DIRS "${CONAN_INCLUDE_DIRS_TCL}")
  set(TCL_INCLUDE_PATH "${CONAN_INCLUDE_DIRS_TCL}")
  set(TCL_LIBRARIES    "${CONAN_PKG_LIBS_TCL}")
  
  message("Win: TCL_LIBRARY:     ${TCL_LIBRARY}")
  message("Win: Env:TCL_LIBRARY: $ENV{TCL_LIBRARY}")
  message("Win: TCL_INCL_PATH:   ${TCL_INCLUDE_PATH}")
  message("Win: TCL:             ${TCL_LIBRARIES}")
# find_package(TCL)
endif()

# pseudo target modeling all exteral packages
add_library(OPS_External_packages INTERFACE)

# include user config
include(${CMAKE_CURRENT_LIST_DIR}/Conf.cmake)

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
    message(":: MacOS")
    add_compile_definitions(_LINUX _UNIX _TCL85)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmax-errors=1")
  else()
    message(":: LINUX")
    add_compile_definitions(_LINUX _UNIX USE_TCL_STUBS _TCL85)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmax-errors=1 -ffloat-store")
    if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address,leak")
    endif()
    add_compile_options(
      $<$<CONFIG:RELEASE>:-march=haswell>
    )
  endif()

  # Debug / Optimization options
  add_compile_options(
    # -pg # Profiling with gprof
    $<$<CONFIG:DEBUG>:-ggdb>
    $<$<CONFIG:DEBUG>:-g3>
    $<$<CONFIG:DEBUG>:-Og>

    $<$<CONFIG:RELEASE>:-O3>
    $<$<CONFIG:RELEASE>:-mtune=native>
    $<$<CONFIG:RELEASE>:-msse4>
  )
  if (DEFINED ProfileBuild)
    add_compile_options(-fno-omit-frame-pointer -g -p)
  endif()
endif()

if(WIN32) # NOTE: this will execute for both 32-bit and 64-bit builds.
  include(OpenSeesDependenciesWin)
  add_compile_definitions(_WIN32 _TCL85)
  add_compile_options(/favor:INTEL64)
  message(STATUS ":: WIN32")
endif()

include(CheckIPOSupported)
check_ipo_supported(RESULT has_ipo OUTPUT error)
if (has_ipo)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()
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
  "$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall>"
  "$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-unused-parameter>"
  "$<$<COMPILE_LANG_AND_ID:CXX,AppleClang,Clang>:-Wno-reorder>"
  "$<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wno-reorder>"
  # $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wsign-conversion>
  # $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wconversion>
)

#----------------------------------------------------------------
# Global Includes
#---------------------------------------------------------------- 

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
  ${OPS_SRC_DIR}/analysis/numberer
  ${OPS_SRC_DIR}/analysis/integrator
  ${OPS_SRC_DIR}/analysis/model

  ${OPS_SRC_DIR}/runtime

  ${OPS_SRC_DIR}/api
  ${OPS_SRC_DIR}/convergenceTest
  ${OPS_SRC_DIR}/domain/constraints

  ${OPS_SRC_DIR}/domain/component
  ${OPS_SRC_DIR}/domain/domain
  ${OPS_SRC_DIR}/database

  ${OPS_SRC_DIR}/domain/domain/single
  ${OPS_SRC_DIR}/domain/groundMotion
  ${OPS_SRC_DIR}/domain/load
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
  ${OPS_SRC_DIR}/recorder
  ${OPS_SRC_DIR}/recorder/response
  ${OPS_SRC_DIR}/system_of_eqn
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE
  ${OPS_SRC_DIR}/system_of_eqn/eigenSOE
  ${OPS_SRC_DIR}/tagged
  ${OPS_SRC_DIR}/tagged/storage
  ${OPS_SRC_DIR}/utility
)

#==============================================================================
#                            Define Targets
#
#==============================================================================
# find_package(MKL)
set(BLA_STATIC ON)
if (TRUE)
  add_subdirectory(OTHER)
else()
  opensees_load(SUPERLU
    #BUNDLED "${OPS_BUNDLED_DIR}/SuperLU_5.1.1/"
    LIBRARY "${OPS_BUNDLED_DIR}/bin/SuperLU_5.1.1/libSUPERLU.a"
  )

  opensees_load(UMFPACK
    #BUNDLED "${OPS_BUNDLED_DIR}/UMFPACK/" 
    LIBRARY "${OPS_BUNDLED_DIR}/bin/UMFPACK/libUMFPACK.a" 
  )

  opensees_load(CSPARSE
    #BUNDLED "${OPS_BUNDLED_DIR}/CSPARSE/"
    LIBRARY "${OPS_BUNDLED_DIR}/bin/CSPARSE/libCSPARSE.a"
  )

  opensees_load(ARPACK
    BUNDLED ${OPS_BUNDLED_DIR}/ARPACK/
  )

  opensees_load(AMD
    BUNDLED "${OPS_BUNDLED_DIR}/AMD/"
    #LIBRARY "${OPS_BUNDLED_DIR}/bin/AMD/libAMD.a"
  )

  opensees_load(METIS
    LIBRARY "${OPS_BUNDLED_DIR}/bin/METIS/libMETIS.a"
  )
endif()

find_package(TCL)
find_package(TclStub) # TODO: may not need to find TCL first
set(TCL_INCLUDE_PATH ${TCL_INCLUDE_DIRS})

include_directories(${TCL_INCLUDE_PATH})

message(":: OpenSees Summary")
message("  :: BLAS:      ${BLAS_LIBRARIES}")
message("  :: LAPACK:    ${LAPACK_LIBRARIES}")
message("  :: SUPERLU:   ${SUPERLU_LIBRARIES}")
message("  :: ARPACK:    ${ARPACK_LIBRARIES}")
message("  :: UMFPACK:   ${UMFPACK_LIBRARIES}")
message("  :: CSPARSE:   ${CSPARSE_LIBRARIES}")
message("  :: AMD:       ${AMD_LIBRARIES}")
message("  :: TCL:       ${TCL_LIBRARY}  ")
message("                ${TCL_STUB_LIBRARY}")
message("                ${TCL_LIBRARIES}")
message("")
message("  :: IPO        ${CMAKE_INTERPROCEDURAL_OPTIMIZATION}")
message("  :: Arch       ${CMAKE_SYSTEM_PROCESSOR}")

# TODO: Probably not needed, link libs as needed
# in sys_of_eqn/linearSOE/*/CMakeLists.txt
add_library(OPS_Numerics           INTERFACE)
target_link_libraries(OPS_Numerics INTERFACE
  ${ARPACK_LIBRARIES} 
  ${SUPERLU_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
)

# Core OpenSees
add_library(OPS_Matrix             OBJECT)
add_library(OPS_Actor              OBJECT)
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
add_library(OPS_Runtime            OBJECT)
add_library(OPS_Database           OBJECT EXCLUDE_FROM_ALL)

# Optional Extensions
add_library(OPS_Parallel           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Parallel_SP        OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_ASDEA              OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Paraview           OBJECT EXCLUDE_FROM_ALL)

# Packaged libraries
add_library(G3)

set_property(TARGET OPS_Runtime           PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET OPS_MaterialFortran   PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET OPS_ElementFortran    PROPERTY POSITION_INDEPENDENT_CODE 1)
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(OPS_ElementFortran PRIVATE "-w")
endif()
set_property(TARGET G3                    PROPERTY POSITION_INDEPENDENT_CODE 1)


# Executables and final targets
#------------------------------------------------------------------------------
# Adding EXCLUDE_FROM_ALL prevents these from being made by default.


# Add sources to targets
add_subdirectory(${OPS_SRC_DIR})

#==============================================================================
#                            Configure targets
#
#==============================================================================
target_link_libraries(OPS_Element      PRIVATE ${OPS_Element_List})
target_link_libraries(OPS_Material     PRIVATE ${TCL_STUB_LIBRARY})


# G3
#----------------------------
target_link_libraries(G3 PUBLIC
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
    OPS_Tagged
    OPS_Utilities
    graph
    coordTransformation
    G3_ObjectBroker
    OPS_Numerics 
    ${OPS_Element_List}
    OPS_Material
    OPS_Damage
    OPS_ModelBuilder
    G3_ObjectBroker
)


#==============================================================================
#                            Apply Options
#
#==============================================================================



# Extensions
#----------------------------
message(":: Configuring OpenSees extensions")
foreach(extension IN LISTS OPS_Element_List OPS_Extension_List)
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

get_target_property(OPS_Damage_COMPILE_OPTIONS OPS_Damage COMPILE_OPTIONS)
  string(REPLACE "-Wall" "" OPS_Damage_COMPILE_OPTIONS "${OPS_Damage_COMPILE_OPTIONS}")
  string(REPLACE "-Wextra" "" OPS_Damage_COMPILE_OPTIONS "${OPS_Damage_COMPILE_OPTIONS}")
  set_target_properties(OPS_Damage PROPERTIES COMPILE_OPTIONS "${OPS_Damage_COMPILE_OPTIONS}")

get_target_property(OPS_Material_COMPILE_OPTIONS OPS_Material COMPILE_OPTIONS)
  string(REPLACE "-Wall" "" OPS_Material_COMPILE_OPTIONS "${OPS_Material_COMPILE_OPTIONS}")
  string(REPLACE "-Wextra" "" OPS_Material_COMPILE_OPTIONS "${OPS_Material_COMPILE_OPTIONS}")
  set_target_properties(OPS_Material PROPERTIES COMPILE_OPTIONS "${OPS_Material_COMPILE_OPTIONS}")

get_target_property(OPS_Element_COMPILE_OPTIONS OPS_Element COMPILE_OPTIONS)
  string(REPLACE "-Wall" "" OPS_Element_COMPILE_OPTIONS "${OPS_Element_COMPILE_OPTIONS}")
  string(REPLACE "-Wextra" "" OPS_Element_COMPILE_OPTIONS "${OPS_Element_COMPILE_OPTIONS}")
  set_target_properties(OPS_Element PROPERTIES COMPILE_OPTIONS "${OPS_Element_COMPILE_OPTIONS}")

