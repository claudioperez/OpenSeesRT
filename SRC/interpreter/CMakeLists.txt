#* ****************************************************************** *#
#*    OpenSees - Open System for Earthquake Engineering Simulation    *#
#*          Pacific Earthquake Engineering Research Center            *#
#* ****************************************************************** *#
cmake_minimum_required(VERSION 3.16)

add_library(OpenSeesRT SHARED)


target_include_directories(OpenSeesRT PUBLIC 
  ./api/
  ${OPS_SRC_DIR}/interpreter/commands
  ${OPS_SRC_DIR}/interpreter/commands/modeling
  ${OPS_SRC_DIR}/interpreter/runtime/modelbuilder
  ${OPS_SRC_DIR}/interpreter/runtime/modelbuilder/basic
  ${OPS_SRC_DIR}/interpreter/runtime/
  ${OPS_SRC_DIR}/interpreter/runtime/Logging/
)

add_subdirectory(commands/)
add_subdirectory(runtime/)

target_sources(OpenSeesRT PRIVATE
# General
    "OpenSeesRT.cpp"
    "G3_Runtime.cpp"
    "elementAPI_PYG3.cpp"
# Other
    "${CMAKE_CURRENT_LIST_DIR}/contrib/packages/optimization/TclParameterCommands.cpp"
#
# SUB-ELEMENT MODELS
#
# Uniaxial
    "${CMAKE_CURRENT_LIST_DIR}/contrib/uniaxial/KikuchiAiken/G3_KikuchiAiken.cpp"
    # "contrib/uniaxial/TclAxialSp.cpp"
    # "contrib/uniaxial/fedeas/TclFedeasWrappers.cpp"
    # "contrib/uniaxial/fedeas/TclFedeasMaterialCommand.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/contrib/uniaxial/PyTzQz/TclPyTzQzMaterialCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/contrib/uniaxial/snap/TclSnapMaterialCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/contrib/uniaxial/drain/TclDrainMaterialCommand.cpp"

    # Other
    #"contrib/uniaxial/fedeas/TclFedeasMaterialCommand.cpp"

#
# OTHER
#
    "${CMAKE_CURRENT_LIST_DIR}/utilities/tclutils.cpp"

)

if (DEFINED OPENSEESRT_VERSION)
  set_property(
      SOURCE OpenSeesRT.cpp
      APPEND PROPERTY COMPILE_DEFINITIONS "OPENSEESRT_VERSION=\"${OPENSEESRT_VERSION}\""
  )
endif()

add_subdirectory(contrib/)

target_compile_definitions(OpenSeesRT PUBLIC USE_TCL_STUBS _TCL85 _NOGRAPHICS)

set_property(TARGET OPS_MaterialFortran            PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET OPS_ElementFortran             PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET OPS_Material_Uniaxial_Drain_f  PROPERTY POSITION_INDEPENDENT_CODE 1)

set_property(TARGET G3         PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET OpenSeesRT PROPERTY POSITION_INDEPENDENT_CODE 1)

if (${TCL_INCLUDE_PATH})
  target_include_directories(OpenSeesRT PUBLIC ${TCL_INCLUDE_PATH})
  target_include_directories(OpenSeesRT_LinearSystem PUBLIC ${TCL_INCLUDE_PATH})
endif()


#
# Python
#
if (NOT DEFINED NoOpenSeesPyRT)
  add_subdirectory(python)
endif()


target_link_libraries(OpenSeesRT PUBLIC ${TCL_LIBRARY})

target_link_libraries(OpenSeesRT PRIVATE
  ${OPS_Element_List}
  OPS_Material
  OPS_Damage
  OPS_ModelBuilder
  G3_ObjectBroker
# OpenSeesRT_LinearSystem
  G3
  ${TCL_STUB_LIBRARY}
  ${UMFPACK_LIBRARIES} ${CSPARSE_LIBRARIES} ${AMD_LIBRARIES}
)

