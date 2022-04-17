#=====================================================================
# 
#   OpenSees -- Open System For Earthquake Engineering Simulation
#           Pacific Earthquake Engineering Research Center
#
#=====================================================================
cmake_minimum_required(VERSION 3.16)
#
# Dependencies
# 
find_package(TclStub)
find_package(TCL)

add_library(OpenSeesRT SHARED)

target_sources(OpenSeesRT PRIVATE
#
# GENERAL COMMANDS
#
# Model Commands & Builders
    "modelbuilder/TclBuilder.cpp"
    "modelbuilder/basic/TclBasicBuilder.cpp"
    #"modelbuilder/plane/TclPlaneStressMaterialTester.cpp"
    "modelbuilder/uniax/TclUniaxialMaterialTester.cpp"
    "modelbuilder/safe/TclSafeBuilder.cpp"
    "modelbuilder/sect/TclSectionTestBuilder.cpp"
    "commands/modeling/geomTransf.cpp"
    "commands/modeling/rigidLink.cpp"

# Other Commands
    "commands/commands.cpp"
    "commands/modeling/model.cpp"
    "commands/modeling/nodes.cpp"
    "commands/analysis/algorithm.cpp"
    "commands/analysis/integrator.cpp"
    "commands/analysis/analysis.cpp"
    "commands/analysis/solver.cpp"
    OpenSeesRT.cpp
    misc/TclPackageClassBroker.cpp

# Database
    "commands/database/TclDatabaseCommands.cpp"
    "commands/database/TclMySQL.cpp"
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

# Multi-Dimensional Materials
    "contrib/nDMaterial/soil/TclUpdateMaterialStageCommand.cpp"
    #"contrib/nD/feap/TclFeapMaterialCommand.cpp"

# Uniaxial
    "contrib/uniaxial/TclKikuchiAiken.cpp"
    "contrib/uniaxial/TclAxialSp.cpp"
    "contrib/uniaxial/fedeas/TclFedeasWrappers.cpp"
    "contrib/uniaxial/fedeas/TclFedeasMaterialCommand.cpp"
    "contrib/uniaxial/PY/TclPyTzQzMaterialCommand.cpp"
    "contrib/uniaxial/limitState/TclLimitState.cpp"
    "contrib/uniaxial/snap/TclSnapMaterialCommand.cpp"
    "contrib/uniaxial/stiffness/TclModelBuilderStiffnessDegradationCommand.cpp"
    "contrib/uniaxial/strength/TclModelBuilderStrengthDegradationCommand.cpp"
    "contrib/uniaxial/unloading/TclModelBuilderUnloadingRuleCommand.cpp"
    "contrib/uniaxial/drain/TclDrainMaterialCommand.cpp"
    "contrib/yieldSurface/evolution/TclModelBuilderYS_EvolutionCommand.cpp"
    "contrib/yieldSurface/plasticHardeningMaterial/TclModelBuilderYSPlasticMaterialCommand.cpp"
    "contrib/yieldSurface/yieldSurfaceBC/TclModelBuilderYieldSurfaceBCCommand.cpp"
    #"material/uniaxial/backbone/TclModelBuilderBackboneCommand.cpp"
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
target_compile_options(OpenSeesRT         PUBLIC  -fPIC)


if (${TCL_INCLUDE_PATH})
  target_include_directories(OpenSeesRT PUBLIC ${TCL_INCLUDE_PATH})
endif()

target_link_libraries(OpenSeesRT PRIVATE
  #${OPS_Extension_List}
  ${OPS_Element_List}
  ${UMFPACK_LIBRARIES} ${CSPARSE_LIBRARIES}
  OPS_Damage
  #OPS_Paraview
  OPS_ModelBuilder
  OPS_Uniaxial_Fedeas
  #OPS_Uniaxial_Fedeas_f
  /home/claudio/packages/FEDEAS/build/src/uniaxial/libfedeas_uniaxial.a
  #/home/claudio/pkgs/FEDEAS/bin/src/uniaxial/libfedeas_uniaxial.a
  G3_ObjectBroker
  #g3_uniaxial_io
  G3
  ${MYSQL_LIBRARIES}
  ${TCL_STUB_LIBRARY}
)


