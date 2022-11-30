#=====================================================================
# 
#   OpenSees -- Open System For Earthquake Engineering Simulation
#           Pacific Earthquake Engineering Research Center
#
#=====================================================================
cmake_minimum_required(VERSION 3.16)

add_library(OpenSeesRT SHARED)

add_library(OpenSeesRT_LinearSystem OBJECT)

target_include_directories(OpenSeesRT PUBLIC 
  ./api/
  ${OPS_SRC_DIR}/api/tclCommandPackage/commands
  ${OPS_SRC_DIR}/api/tclCommandPackage/commands/modeling
  ${OPS_SRC_DIR}/api/tclCommandPackage/misc
  ${OPS_SRC_DIR}/api/tclCommandPackage/modelbuilder
  ${OPS_SRC_DIR}/api/tclCommandPackage/modelbuilder/basic
  ${OPS_SRC_DIR}/api/tclCommandPackage/modelbuilder/uniax
  ${OPS_SRC_DIR}/api/tclCommandPackage/modelbuilder/sect
  ${OPS_SRC_DIR}/api/tclCommandPackage/runtime/
  ${OPS_SRC_DIR}/api/tclCommandPackage/runtime/Logging/
)

add_subdirectory(commands/)
add_subdirectory(runtime/)

target_sources(OpenSeesRT PRIVATE
#
# GENERAL
#
    "OpenSeesRT.cpp"
    "G3_Runtime.cpp"
    "elementAPI_PYG3.cpp"

# Model Builders
    "modelbuilder/TclBuilder.cpp"
    "modelbuilder/basic/TclBasicBuilder.cpp"
    "modelbuilder/uniax/TclUniaxialMaterialTester.cpp"
    "modelbuilder/sect/TclSectionTestBuilder.cpp"
    #"modelbuilder/plane/TclPlaneStressMaterialTester.cpp"


# Other
    "misc/TclPackageClassBroker.cpp"
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
    #"contrib/uniaxial/limitState/TclLimitState.cpp"
    "contrib/uniaxial/snap/TclSnapMaterialCommand.cpp"
    "contrib/uniaxial/drain/TclDrainMaterialCommand.cpp"

    # Other
    #"contrib/uniaxial/fedeas/TclFedeasMaterialCommand.cpp"

#
# OTHER
#
    # "ffi/C_Wrapper.cpp"
    "misc/tcl_packages.cpp"
    "utilities/tclutils.cpp"

)

add_subdirectory(contrib/)

target_compile_definitions(OpenSeesRT PUBLIC USE_TCL_STUBS _TCL85 _NOGRAPHICS)

target_compile_options(OPS_MaterialFortran            PRIVATE -fPIC)
target_compile_options(OPS_ElementFortran             PRIVATE -fPIC)
target_compile_options(OPS_Material_Uniaxial_Drain_f  PRIVATE -fPIC)
set_property(TARGET G3         PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET OpenSeesRT PROPERTY POSITION_INDEPENDENT_CODE 1)

if (${TCL_INCLUDE_PATH})
  target_include_directories(OpenSeesRT PUBLIC ${TCL_INCLUDE_PATH})
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
  ${UMFPACK_LIBRARIES} ${CSPARSE_LIBRARIES} ${AMD_LIBRARIES}
  OPS_Material
  OPS_Damage
  OPS_ModelBuilder
  G3_ObjectBroker
  OpenSeesRT_LinearSystem
  G3
  ${TCL_STUB_LIBRARY}
)

