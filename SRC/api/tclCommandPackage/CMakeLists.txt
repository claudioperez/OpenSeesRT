#=====================================================================
# 
#   OpenSees -- Open System For Earthquake Engineering Simulation
#           Pacific Earthquake Engineering Research Center
#
#=====================================================================
cmake_minimum_required(VERSION 3.16)

add_library(OpenSeesRT SHARED)

# add_library(OpenSeesRT_Integrator   OBJECT)
# target_sources(OpenSeesRT_Integrator PRIVATE
#     "commands/analysis/integrator.cpp"
# )

add_library(OpenSeesRT_LinearSystem OBJECT)
target_sources(OpenSeesRT_LinearSystem PRIVATE
    "commands/analysis/solver.cpp"
    "commands/analysis/solver.hpp"
)

target_sources(OpenSeesRT PRIVATE
#
# GENERAL COMMANDS
#
# Model Commands & Builders
    "G3_Runtime.cpp"
    "modelbuilder/TclBuilder.cpp"
    "modelbuilder/basic/TclBasicBuilder.cpp"
    #"modelbuilder/plane/TclPlaneStressMaterialTester.cpp"
    "modelbuilder/uniax/TclUniaxialMaterialTester.cpp"
    "modelbuilder/safe/TclSafeBuilder.cpp"
    "modelbuilder/sect/TclSectionTestBuilder.cpp"
    "commands/modeling/geomTransf.cpp"
    "commands/modeling/rigidLink.cpp"
    "commands/domain.cpp"

# Other Commands
    "streams/logging.cpp"
    "commands/commands.cpp"
    "commands/runtime.cpp"
    "commands/modeling/model.cpp"
    "commands/modeling/nodes.cpp"
    "commands/modeling/constraint.cpp"
    "commands/analysis/algorithm.cpp"
    "commands/analysis/integrator.cpp"
    "commands/analysis/analysis.cpp"
    "commands/analysis/ctest.cpp"

    OpenSeesRT.cpp
    misc/TclPackageClassBroker.cpp
    "commands/recorder.cpp"
    "contrib/packages/optimization/TclParameterCommands.cpp"

# Database
    # "commands/database/TclDatabaseCommands.cpp"
    # "commands/database/TclMySQL.cpp"
    #"commands/database/TclBerkeleyDB.cpp"

# Multi-Dimensional Materials
    #"contrib/nD/feap/TclFeapMaterialCommand.cpp"
    "commands/modeling/TclUpdateMaterialCommand.cpp"
    "commands/TclUpdateMaterialStageCommand.cpp"

# LOADS & PATTERNS
    "commands/loading/groundMotion.cpp"
    "commands/loading/TclPatternCommand.cpp"
    "commands/loading/TclSeriesCommand.cpp"
    "commands/loading/TclSeriesIntegratorCommand.cpp"
    #"domain/pattern/drm/TclPatternCommand.cpp"
    "commands/modeling/region.cpp"
#
# ELEMENTS
#
    "commands/modeling/element.cpp"
    #"element/feap/TclFeapElementCommand.cpp"
    #"element/PFEMElement/TclModelBuilder_addPFEMElement.cpp"
    #"element/PFEMElement/TclPFEMCommands.cpp"

#
# SUB-ELEMENT MODELS
#
# Commands
    "commands/modeling/nDMaterial.cpp"
    "commands/modeling/section.cpp"
    "commands/modeling/uniaxialMaterial.cpp"
    "commands/modeling/uniaxial.cpp"


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
    "contrib/uniaxial/limitState/TclLimitState.cpp"
    "contrib/uniaxial/snap/TclSnapMaterialCommand.cpp"
    "contrib/uniaxial/drain/TclDrainMaterialCommand.cpp"

    # Other
    #"material/uniaxial/TclNewUniaxialMaterial.cpp"
    #"contrib/uniaxial/fedeas/TclFedeasMaterialCommand.cpp"


#
# OTHER
#
    #"contrib/system_of_eqn/TclPetsc.cpp"

    "elementAPI_PYG3.cpp"
    "api/C_Wrapper.cpp"
    "misc/tcl_packages.cpp"
    "misc/TclFeViewer.cpp"
    "utilities/tclutils.cpp"
)

# add_subdirectory(libraries/uniaxial)
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
  # target_link_libraries(OpenSeesRT PUBLIC OpenSeesPyRT)
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

