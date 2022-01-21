#=====================================================================
# 
#   OpenSees -- Open System For Earthquake Engineering Simulation
#           Pacific Earthquake Engineering Research Center
#
#=====================================================================
cmake_minimum_required(VERSION 3.16)

# Dependencies
find_package(TclStub)
find_package(TCL)

add_library(OpenSeesCommandPackage SHARED)

target_sources(OpenSeesCommandPackage PRIVATE
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
    "commands/myCommands.cpp"

    "commands/modeling/nodes.cpp"
    "commands/analysis/algorithm.cpp"
    "commands/analysis/integrator.cpp"
    "commands/analysis/analysis.cpp"
    OpenSeesCommandPackage.cpp
    TclPackageClassBroker.cpp

# Database
    "database/TclDatabaseCommands.cpp"
    "database/TclMySQL.cpp"
    #"database/TclBerkeleyDB.cpp"

    "recorder/TclRecorderCommands.cpp"
    "domain/component/TclParameterCommands.cpp"
    "domain/component/TclUpdateMaterialCommand.cpp"
# LOADS & PATTERNS
    "domain/groundMotion/TclGroundMotionCommand.cpp"
    "domain/pattern/TclPatternCommand.cpp"
    "domain/pattern/TclSeriesCommand.cpp"
    "domain/pattern/TclSeriesIntegratorCommand.cpp"
    #"domain/pattern/drm/TclPatternCommand.cpp"
    "domain/region/TclRegionCommands.cpp"
#
# ELEMENTS
#
    "commands/modeling/element.cpp"
    "element/beamWithHinges.cpp"
    "element/gradientInelasticBeamColumn.cpp"
    "element/UP-ucsd/TclFourNodeQuadUPCommand.cpp"
    #"element/UWelements/Tcl_generateInterfacePoints.cpp"
    "element/adapter/TclActuatorCommand.cpp"
    "element/adapter/TclActuatorCorotCommand.cpp"
    "element/adapter/TclAdapterCommand.cpp"
    "element/brick/TclBrickCommand.cpp"
    "element/brick/TclTwenty_Node_BrickCommand.cpp"
    "element/dispBeamColumnInt/TclDispBeamColumnIntCommand.cpp"
    "element/elasticBeamColumn/TclElasticBeamCommand.cpp"
    "element/elastomericBearing/TclElastomericBearingBoucWenCommand.cpp"
    "element/elastomericBearing/TclElastomericBearingPlasticityCommand.cpp"
    "element/elastomericBearing/TclElastomericBearingUFRPCommand.cpp"
    "element/forceBeamColumn/TclForceBeamColumnCommand.cpp"
    "element/fourNodeQuad/TclFourNodeQuadCommand.cpp"
    "element/frictionBearing/TclFlatSliderCommand.cpp"
    "element/frictionBearing/TclRJWatsonEQSCommand.cpp"
    "element/frictionBearing/TclSingleFPCommand.cpp"
    "element/frictionBearing/frictionModel/TclModelBuilderFrictionModelCommand.cpp"
    "element/generic/TclGenericClientCommand.cpp"
    "element/generic/TclGenericCopyCommand.cpp"
    "element/joint/TclBeamColumnJointCommand.cpp"
    "element/joint/TclJoint2dCommand.cpp"
    "element/joint/TclJoint3dCommand.cpp"
    "element/twoNodeLink/TclTwoNodeLinkCommand.cpp"
    "element/updatedLagrangianBeamColumn/TclCyclicModelCommands.cpp"
    "element/updatedLagrangianBeamColumn/TclElement2dGNL.cpp"
    "element/updatedLagrangianBeamColumn/TclElement2dYS.cpp"
    "element/zeroLength/TclZeroLength.cpp"
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
    "material/nD/soil/TclUpdateMaterialStageCommand.cpp"
    #"material/nD/feap/TclFeapMaterialCommand.cpp"

# Uniaxial
    "contrib/uniaxial/TclKikuchiAiken.cpp"
    "contrib/uniaxial/TclAxialSp.cpp"
    "contrib/uniaxial/fedeas/TclFedeasWrappers.cpp"
    "material/uniaxial/PY/TclPyTzQzMaterialCommand.cpp"
    "material/uniaxial/limitState/TclLimitState.cpp"
    "material/uniaxial/snap/TclSnapMaterialCommand.cpp"
    "material/uniaxial/stiffness/TclModelBuilderStiffnessDegradationCommand.cpp"
    "material/uniaxial/strength/TclModelBuilderStrengthDegradationCommand.cpp"
    "material/uniaxial/unloading/TclModelBuilderUnloadingRuleCommand.cpp"
    "material/uniaxial/drain/TclDrainMaterialCommand.cpp"
    "material/yieldSurface/evolution/TclModelBuilderYS_EvolutionCommand.cpp"
    "material/yieldSurface/plasticHardeningMaterial/TclModelBuilderYSPlasticMaterialCommand.cpp"
    "material/yieldSurface/yieldSurfaceBC/TclModelBuilderYieldSurfaceBCCommand.cpp"
    #"material/uniaxial/backbone/TclModelBuilderBackboneCommand.cpp"
    #"material/uniaxial/TclNewUniaxialMaterial.cpp"
    #"contrib/uniaxial/fedeas/TclFedeasMaterialCommand.cpp"

# Section
    "material/section/yieldSurface/TclModelBuilderYS_SectionCommand.cpp"

#
# Reliability & Optimization 
#
    #"optimization/tcl/TclOptimizationBuilder.cpp"

    #"reliability/analysis/gFunction/TclGFunEvaluator.cpp"
    #"reliability/analysis/gFunction/TclMatlabGFunEvaluator.cpp"
    #"reliability/domain/functionEvaluator/TclEvaluator.cpp"
    #"reliability/tcl/TclReliabilityBuilder.cpp"
#
# OTHER
#
    #"system_of_eqn/linearSOE/petsc/TclPetsc.cpp"
    
    #"${OPS_SRC_DIR}/interpreter/OpenSeesMiscCommands.cpp"

    "elementAPI_PYG3.cpp"
    "packages/tcl_packages.cpp"
    "misc/TclFeViewer.cpp"
)

add_subdirectory(libraries/uniaxial)

target_compile_definitions(OpenSeesCommandPackage PUBLIC USE_TCL_STUBS _TCL85 _NOGRAPHICS)

target_compile_options(OPS_MaterialFortran            PRIVATE -fPIC)
target_compile_options(OPS_ElementFortran             PRIVATE -fPIC)
target_compile_options(OPS_Material_Uniaxial_Drain_f  PRIVATE -fPIC)
target_compile_options(OpenSeesCommandPackage         PUBLIC  -fPIC)


if (${TCL_INCLUDE_PATH})
target_include_directories(OpenSeesCommandPackage PUBLIC ${TCL_INCLUDE_PATH})
endif()

target_link_libraries(OpenSeesCommandPackage PRIVATE
  #${OPS_Extension_List}
  ${OPS_Element_List}
  ${UMFPACK_LIBRARIES} ${CSPARSE_LIBRARIES}
  OPS_Damage
  #OPS_Paraview
  OPS_ModelBuilder
  G3_ObjectBroker
  #g3_uniaxial_io
  G3
  ${MYSQL_LIBRARIES}
  ${TCL_STUB_LIBRARY}
)


