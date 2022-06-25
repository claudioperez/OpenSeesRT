
add_subdirectory(uniaxial)
#add_subdirectory(system_of_eqn)
# add_subdirectory(packages/YieldSurface/)

target_sources(OpenSeesRT PRIVATE
    "integrator/Newmark/G3Parse_Newmark.cpp"

    "element/UP-ucsd/TclFourNodeQuadUPCommand.cpp"
    #"element/UWelements/Tcl_generateInterfacePoints.cpp"
    "element/adapter/TclActuatorCommand.cpp"
    "element/adapter/TclActuatorCorotCommand.cpp"
    "element/adapter/TclAdapterCommand.cpp"

    "element/brick/TclBrickCommand.cpp"
    "element/brick/TclTwenty_Node_BrickCommand.cpp"
    "element/elasticBeamColumn/TclElasticBeamCommand.cpp"
    "element/elastomericBearing/TclElastomericBearingBoucWenCommand.cpp"
    "element/elastomericBearing/TclElastomericBearingPlasticityCommand.cpp"
    "element/elastomericBearing/TclElastomericBearingUFRPCommand.cpp"
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
    "element/zeroLength/TclZeroLength.cpp"

    "element/nonlinearBeamColumn/updatedLagrangianBeamColumn/TclCyclicModelCommands.cpp"
    "element/nonlinearBeamColumn/updatedLagrangianBeamColumn/TclElement2dGNL.cpp"
    "element/nonlinearBeamColumn/updatedLagrangianBeamColumn/TclElement2dYS.cpp"
    "element/nonlinearBeamColumn/forceBeamColumn/TclForceBeamColumnCommand.cpp"
    "element/nonlinearBeamColumn/dispBeamColumnInt/TclDispBeamColumnIntCommand.cpp"
    "element/nonlinearBeamColumn/beamWithHinges.cpp"
    "element/nonlinearBeamColumn/gradientInelasticBeamColumn.cpp"
)


