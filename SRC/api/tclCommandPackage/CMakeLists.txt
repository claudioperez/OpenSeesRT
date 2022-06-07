#=====================================================================
# 
#   OpenSees -- Open System For Earthquake Engineering Simulation
#           Pacific Earthquake Engineering Research Center
#
#=====================================================================
cmake_minimum_required(VERSION 3.16)

add_library(OpenSeesRT SHARED)

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
    "commands/analysis/solver.cpp"
    OpenSeesRT.cpp
    misc/TclPackageClassBroker.cpp

# Database
# "commands/database/TclDatabaseCommands.cpp"
# "commands/database/TclMySQL.cpp"
    #"commands/database/TclBerkeleyDB.cpp"

    "commands/recorder.cpp"
    "commands/modeling/TclParameterCommands.cpp"
    "commands/modeling/TclUpdateMaterialCommand.cpp"

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
    #"element/nonlinearBeamColumn/tcl/TclElmtBuilder.cpp"
#
# SUB-ELEMENT MODELS
#
# Commands
    "commands/modeling/nDMaterial.cpp"
    "commands/modeling/section.cpp"
    "commands/modeling/uniaxialMaterial.cpp"
    "commands/modeling/uniaxial.cpp"
    "commands/modeling/hystereticBackbone.cpp"
    "commands/modeling/backbone.cpp"

# Multi-Dimensional Materials
    "contrib/nDMaterial/soil/TclUpdateMaterialStageCommand.cpp"
    #"contrib/nD/feap/TclFeapMaterialCommand.cpp"

# Uniaxial
    "contrib/uniaxial/KikuchiAiken/G3_KikuchiAiken.cpp"
    # "contrib/uniaxial/TclAxialSp.cpp"
    # "contrib/uniaxial/fedeas/TclFedeasWrappers.cpp"
    # "contrib/uniaxial/fedeas/TclFedeasMaterialCommand.cpp"
    # "contrib/uniaxial/ConcreteXX/concrete04.cpp"
    # "contrib/uniaxial/ReinforcingSteel/TclReinforcingSteel.cpp"

    "contrib/uniaxial/PyTzQz/TclPyTzQzMaterialCommand.cpp"
    "contrib/uniaxial/limitState/TclLimitState.cpp"
    "contrib/uniaxial/snap/TclSnapMaterialCommand.cpp"
    "contrib/uniaxial/drain/TclDrainMaterialCommand.cpp"
    # YieldSurface
    "contrib/yieldSurface/TclModelBuilderYS_EvolutionCommand.cpp"
    "contrib/yieldSurface/TclModelBuilderYSPlasticMaterialCommand.cpp"
    "contrib/yieldSurface/TclModelBuilderYieldSurfaceBCCommand.cpp"
    # OO_Hysteretic
    "contrib/uniaxial/OO_Hysteretic/StiffnessDegradationCommand.cpp"
    "contrib/uniaxial/OO_Hysteretic/StrengthDegradationCommand.cpp"
    "contrib/uniaxial/OO_Hysteretic/UnloadingRuleCommand.cpp"

    #"material/uniaxial/backbone/TclModelBuilderBackboneCommand.cpp"
    # Other
    #"material/uniaxial/TclNewUniaxialMaterial.cpp"
    #"contrib/uniaxial/fedeas/TclFedeasMaterialCommand.cpp"

# Section
    "contrib/yieldSurface/TclModelBuilderYS_SectionCommand.cpp"

#
# Reliability & Optimization 
#
    #"contrib/optimization/TclOptimizationBuilder.cpp"

    #"contrib/reliability/analysis/gFunction/TclGFunEvaluator.cpp"
    #"contrib/reliability/analysis/gFunction/TclMatlabGFunEvaluator.cpp"
    #"contrib/reliability/domain/functionEvaluator/TclEvaluator.cpp"
    #"contrib/reliability/tcl/TclReliabilityBuilder.cpp"
#
# OTHER
#
    #"contrib/system_of_eqn/TclPetsc.cpp"
    
    #"${OPS_SRC_DIR}/interpreter/OpenSeesMiscCommands.cpp"

    "elementAPI_PYG3.cpp"
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
# FEDEAS uniaxial materials
#
set(fedeas_uniaxial "/home/claudio/packages/FEDEAS/build/src/uniaxial/libfedeas_uniaxial.a")
if (EXISTS "${fedeas_uniaxial}")
  target_link_libraries(OpenSeesRT PRIVATE
    OPS_Uniaxial_Fedeas
    ${fedeas_uniaxial}
  )
endif()

#
# Python
#
if (NOT DEFINED NoOpenSeesPyRT)
  add_subdirectory(python)
  target_link_libraries(OpenSeesRT PRIVATE OpenSeesPyRT)
endif()


target_link_libraries(OpenSeesRT PRIVATE
  #${OPS_Extension_List}
  ${OPS_Element_List}
  ${UMFPACK_LIBRARIES} ${CSPARSE_LIBRARIES}
  OPS_Damage
  #OPS_Paraview
  OPS_ModelBuilder
  G3_ObjectBroker
  G3
  # OpenSeesPyRT
  # ${MYSQL_LIBRARIES}
  ${TCL_STUB_LIBRARY}
  #${TCL_LIBRARIES}
)

