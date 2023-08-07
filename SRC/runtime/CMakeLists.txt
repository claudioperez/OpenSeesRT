#* ****************************************************************** *#
#*    OpenSees - Open System for Earthquake Engineering Simulation    *#
#*          Pacific Earthquake Engineering Research Center            *#
#* ****************************************************************** *#

add_library(OpenSeesRT SHARED)
target_sources(OpenSeesRT PRIVATE "OpenSeesRT.cpp")
target_compile_definitions(OpenSeesRT PUBLIC USE_TCL_STUBS _TCL85)
set_property(TARGET OpenSeesRT            PROPERTY POSITION_INDEPENDENT_CODE 1)

target_compile_options(OPS_Runtime PRIVATE
  "$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wextra>"
  "$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wpedantic>"
  "$<$<CXX_COMPILER_ID:MSVC>:/W4>"
)

target_include_directories(OPS_Runtime PUBLIC 
  ./api/
  ${OPS_SRC_DIR}/runtime/commands
  ${OPS_SRC_DIR}/runtime/commands/modeling
  ${OPS_SRC_DIR}/runtime/runtime/modelbuilder
  ${OPS_SRC_DIR}/runtime/runtime/modelbuilder/basic
  ${OPS_SRC_DIR}/runtime/runtime/
  ${OPS_SRC_DIR}/runtime/runtime/Logging/
)
target_link_libraries(OPS_Runtime PRIVATE G3)

add_subdirectory(commands/)

add_subdirectory(runtime/)

target_sources(OPS_Runtime PRIVATE
    "G3_Runtime.cpp"
    "elementAPI_PYG3.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/utilities/tclutils.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/contrib/packages/optimization/TclParameterCommands.cpp"

    #
    # OTHER
    #
    "${CMAKE_CURRENT_LIST_DIR}/contrib/uniaxial/KikuchiAiken/G3_KikuchiAiken.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/contrib/uniaxial/PyTzQz/TclPyTzQzMaterialCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/contrib/uniaxial/snap/TclSnapMaterialCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/contrib/uniaxial/drain/TclDrainMaterialCommand.cpp"
    # "contrib/uniaxial/TclAxialSp.cpp"
    # "contrib/uniaxial/fedeas/TclFedeasWrappers.cpp"
    # "contrib/uniaxial/fedeas/TclFedeasMaterialCommand.cpp"
    #"contrib/uniaxial/fedeas/TclFedeasMaterialCommand.cpp"

)

if (DEFINED OPENSEESRT_VERSION)
  set_property(
    SOURCE OpenSeesRT.cpp
    APPEND PROPERTY COMPILE_DEFINITIONS 
           "OPENSEESRT_VERSION=\"${OPENSEESRT_VERSION}\""
  )
endif()

add_subdirectory(contrib/)

#
# Tcl
#
if (${TCL_INCLUDE_PATH})
  target_include_directories(OpenSeesRT  PUBLIC ${TCL_INCLUDE_PATH})
  target_include_directories(OPS_Runtime PUBLIC ${TCL_INCLUDE_PATH})
endif()

target_link_libraries(OpenSeesRT PRIVATE OPS_Runtime ${TCL_STUB_LIBRARY})

#
# Python
#
  message(" :: Configuring Pyton extension")
if (NOT DEFINED NoOpenSeesPyRT)
  add_subdirectory(python)
endif()

#
# Parallel
#
if (MPI_FOUND)
  add_subdirectory(parallel)
endif()

