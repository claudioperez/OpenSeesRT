#=====================================================================
# 
#   OpenSees -- Open System For Earthquake Engineering Simulation
#           Pacific Earthquake Engineering Research Center
#
#=====================================================================
cmake_minimum_required(VERSION 3.16)

add_library(OpenSeesRT SHARED)

add_library(OpenSeesRT_LinearSystem OBJECT)

target_sources(OpenSeesRT_LinearSystem 
  PRIVATE
    "commands/analysis/solver.cpp"
  PUBLIC
    "commands/analysis/solver.hpp"
)

add_subdirectory(commands/)

target_sources(OpenSeesRT PRIVATE
#
# GENERAL COMMANDS
#
# Model Builders
    "OpenSeesRT.cpp"
    "G3_Runtime.cpp"
    "streams/logging.cpp"
    "modelbuilder/TclBuilder.cpp"
    "modelbuilder/basic/TclBasicBuilder.cpp"
    "modelbuilder/uniax/TclUniaxialMaterialTester.cpp"
    "modelbuilder/safe/TclSafeBuilder.cpp"
    "modelbuilder/sect/TclSectionTestBuilder.cpp"
    "runtime/BasicAnalysisBuilder.cpp"
    #"modelbuilder/plane/TclPlaneStressMaterialTester.cpp"


# Other
    "misc/TclPackageClassBroker.cpp"
    "contrib/packages/optimization/TclParameterCommands.cpp"

#
# SUB-ELEMENT MODELS
#
    # YIELD SURFACE
    "contrib/packages/YieldSurface/TclModelBuilderYS_EvolutionCommand.cpp"
    "contrib/packages/YieldSurface/TclModelBuilderYSPlasticMaterialCommand.cpp"
    "contrib/packages/YieldSurface/TclModelBuilderYieldSurfaceBCCommand.cpp"
    "contrib/packages/YieldSurface/TclModelBuilderYS_SectionCommand.cpp"

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
    "elementAPI_PYG3.cpp"
    # "ffi/C_Wrapper.cpp"
    "misc/tcl_packages.cpp"
    "utilities/tclutils.cpp"

)

add_subdirectory(contrib/)

target_compile_definitions(OpenSeesRT PUBLIC USE_TCL_STUBS _TCL85 _NOGRAPHICS)

target_compile_options(OPS_MaterialFortran            PRIVATE -fPIC)
target_compile_options(OPS_ElementFortran             PRIVATE -fPIC)
target_compile_options(OPS_Material_Uniaxial_Drain_f  PRIVATE -fPIC)
target_compile_options(OpenSeesRT                     PUBLIC  -fPIC)


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
  OPS_Damage
  OPS_ModelBuilder
  G3_ObjectBroker
  OpenSeesRT_LinearSystem
  G3
  # ${MYSQL_LIBRARIES}
  ${TCL_STUB_LIBRARY}
)

