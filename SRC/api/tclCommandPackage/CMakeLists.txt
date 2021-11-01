#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
cmake_minimum_required(VERSION 3.16)


add_library(OpenSeesCommandPackage SHARED)

target_sources(OpenSeesCommandPackage
  PRIVATE
    # Model Builders
    "modelbuilder/basic/TclBasicBuilder.cpp"
    #"modelbuilder/plane/TclPlaneStressMaterialTester.cpp"
    #"modelbuilder/uniax/TclUniaxialMaterialTester.cpp"
    "modelbuilder/safe/TclSafeBuilder.cpp"
    #"modelbuilder/sect/TclSectionTester.cpp"

    "commands.cpp"
    "myCommands.cpp"
    OpenSeesCommandPackage.cpp
    TclPackageClassBroker.cpp
    "coordTransformation/TclGeomTransfCommand.cpp"

    ## Database
    "database/TclDatabaseCommands.cpp"
    "database/TclMySQL.cpp"
    #"database/TclBerkeleyDB.cpp"

    "recorder/TclRecorderCommands.cpp"
    "domain/component/TclParameterCommands.cpp"
    "domain/component/TclUpdateMaterialCommand.cpp"
    "domain/groundMotion/TclGroundMotionCommand.cpp"
    "domain/pattern/TclPatternCommand.cpp"
    "domain/pattern/TclSeriesCommand.cpp"
    "domain/pattern/TclSeriesIntegratorCommand.cpp"
    #"domain/pattern/drm/TclPatternCommand.cpp"
    "domain/region/TclRegionCommands.cpp"

    "element/TclElementCommands.cpp"
    "element/UP-ucsd/TclFourNodeQuadUPCommand.cpp"
    "element/UWelements/Tcl_generateInterfacePoints.cpp"
    "element/adapter/TclActuatorCommand.cpp"
    "element/adapter/TclActuatorCorotCommand.cpp"
    "element/adapter/TclAdapterCommand.cpp"
    "element/beamWithHinges/TclBeamWithHingesBuilder.cpp"
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
    "element/gradientInelasticBeamColumn/TclGradientInelasticBeamColumnCommand.cpp"
    "element/joint/TclBeamColumnJointCommand.cpp"
    "element/joint/TclJoint2dCommand.cpp"
    "element/joint/TclJoint3dCommand.cpp"
    "element/twoNodeLink/TclTwoNodeLinkCommand.cpp"
    "element/updatedLagrangianBeamColumn/TclCyclicModelCommands.cpp"
    "element/updatedLagrangianBeamColumn/TclElement2dGNL.cpp"
    "element/updatedLagrangianBeamColumn/TclElement2dYS.cpp"
    "element/zeroLength/TclZeroLength.cpp"
    #"element/nonlinearBeamColumn/tcl/TclElmtBuilder.cpp"

    ## Materials
    "material/nD/TclModelBuilderNDMaterialCommand.cpp"
    "material/nD/soil/TclUpdateMaterialStageCommand.cpp"

    "material/section/TclModelBuilderSectionCommand.cpp"
    "material/section/yieldSurface/TclModelBuilderYS_SectionCommand.cpp"

    "material/uniaxial/PY/TclPyTzQzMaterialCommand.cpp"
    "material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp"
    "material/uniaxial/backbone/TclModelBuilderBackboneCommand.cpp"
    "material/uniaxial/limitState/TclLimitState.cpp"
    "material/uniaxial/snap/TclSnapMaterialCommand.cpp"
    "material/uniaxial/stiffness/TclModelBuilderStiffnessDegradationCommand.cpp"
    "material/uniaxial/strength/TclModelBuilderStrengthDegradationCommand.cpp"
    "material/uniaxial/unloading/TclModelBuilderUnloadingRuleCommand.cpp"
    #"material/uniaxial/TclNewUniaxialMaterial.cpp"
    "material/yieldSurface/evolution/TclModelBuilderYS_EvolutionCommand.cpp"
    "material/yieldSurface/plasticHardeningMaterial/TclModelBuilderYSPlasticMaterialCommand.cpp"
    "material/yieldSurface/yieldSurfaceBC/TclModelBuilderYieldSurfaceBCCommand.cpp"

    #"damage/TclModelBuilderDamageModelCommand.cpp"

    "elementAPI_PYG3.cpp"
    tcl_packages.cpp
    "misc/TclFeViewer.cpp"

    ## Optimization 
    #"optimization/tcl/TclOptimizationBuilder.cpp"

    #"reliability/analysis/gFunction/TclGFunEvaluator.cpp"
    #"reliability/analysis/gFunction/TclMatlabGFunEvaluator.cpp"
    #"reliability/domain/functionEvaluator/TclEvaluator.cpp"
    #"reliability/tcl/TclReliabilityBuilder.cpp"
    #"system_of_eqn/linearSOE/petsc/TclPetsc.cpp"
    #"tcl/TclVideoPlayer.cpp"
    
    #"material/nD/feap/TclFeapMaterialCommand.cpp"
    "material/uniaxial/drain/TclDrainMaterialCommand.cpp"
    #"material/uniaxial/fedeas/TclFedeasMaterialCommand.cpp"
    #"element/feap/TclFeapElementCommand.cpp"
    #"element/PFEMElement/TclModelBuilder_addPFEMElement.cpp"
    #"element/PFEMElement/TclPFEMCommands.cpp"
    
    #"interpreter/TclWrapper.cpp"
    #"${OPS_SRC_DIR}/interpreter/OpenSeesMiscCommands.cpp"
)

target_compile_definitions(OpenSeesCommandPackage PRIVATE USE_TCL_STUBS _TCL85 _NOGRAPHICS)
target_compile_options(OpenSeesCommandPackage PUBLIC -fPIC)

# Dependencies
find_package(TclStub)
find_package(TCL)

if (${TCL_INCLUDE_PATH})
target_include_directories(OpenSeesCommandPackage PUBLIC ${TCL_INCLUDE_PATH})
endif()

target_link_libraries(OpenSeesCommandPackage PRIVATE
  #${OPS_Extension_List}
  ${OPS_Element_List}
  #-Wl,--whole-archive
  #OPS_Uniaxial_Fedeas
  OPS_Element_PFEMElement OPS_PFEM ${UMFPACK_LIBRARIES} ${CSPARSE_LIBRARIES}
  #OpenSees 
  OPS_Damage
  OPS_Paraview
  OPS_ModelBuilder
  G3_ObjectBroker
  G3
  ${MYSQL_LIBRARIES}
  ${TCL_STUB_LIBRARY}
)


