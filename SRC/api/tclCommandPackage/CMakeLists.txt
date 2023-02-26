#* ****************************************************************** *#
#*    OpenSees - Open System for Earthquake Engineering Simulation    *#
#*          Pacific Earthquake Engineering Research Center            *#
#* ****************************************************************** *#
cmake_minimum_required(VERSION 3.16)

add_library(OpenSeesRT SHARED)

add_library(OpenSeesRT_LinearSystem OBJECT)
target_compile_definitions(OpenSeesRT_LinearSystem PUBLIC USE_TCL_STUBS)

target_include_directories(OpenSeesRT PUBLIC 
  ./api/
  ${OPS_SRC_DIR}/api/tclCommandPackage/commands
  ${OPS_SRC_DIR}/api/tclCommandPackage/commands/modeling
  ${OPS_SRC_DIR}/api/tclCommandPackage/runtime/modelbuilder
  ${OPS_SRC_DIR}/api/tclCommandPackage/runtime/modelbuilder/basic
  ${OPS_SRC_DIR}/api/tclCommandPackage/runtime/modelbuilder/uniax
  ${OPS_SRC_DIR}/api/tclCommandPackage/runtime/modelbuilder/sect
  ${OPS_SRC_DIR}/api/tclCommandPackage/runtime/
  ${OPS_SRC_DIR}/api/tclCommandPackage/runtime/Logging/
)

add_subdirectory(commands/)
add_subdirectory(runtime/)

target_sources(OpenSeesRT PRIVATE
# General
    "OpenSeesRT.cpp"
    "G3_Runtime.cpp"
    "elementAPI_PYG3.cpp"

# Other
    "contrib/packages/optimization/TclParameterCommands.cpp"

#
# SUB-ELEMENT MODELS
#
# Uniaxial
    "contrib/uniaxial/KikuchiAiken/G3_KikuchiAiken.cpp"
    # "contrib/uniaxial/TclAxialSp.cpp"
    # "contrib/uniaxial/fedeas/TclFedeasWrappers.cpp"
    # "contrib/uniaxial/fedeas/TclFedeasMaterialCommand.cpp"

    "contrib/uniaxial/PyTzQz/TclPyTzQzMaterialCommand.cpp"
    "contrib/uniaxial/snap/TclSnapMaterialCommand.cpp"
    "contrib/uniaxial/drain/TclDrainMaterialCommand.cpp"

    # Other
    #"contrib/uniaxial/fedeas/TclFedeasMaterialCommand.cpp"

#
# OTHER
#
    "utilities/tclutils.cpp"

)

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

target_compile_options(OpenSeesRT PRIVATE -fmax-errors=1)


#
# Python
#
if (NOT DEFINED NoOpenSeesPyRT)
  add_subdirectory(python)
endif()


target_link_libraries(OpenSeesRT PRIVATE
  ${OPS_Element_List}
  OPS_Material
  OPS_Damage
  OPS_ModelBuilder
  G3_ObjectBroker
  OpenSeesRT_LinearSystem
  G3
  ${TCL_STUB_LIBRARY}
  ${UMFPACK_LIBRARIES} ${CSPARSE_LIBRARIES} ${AMD_LIBRARIES}
)

