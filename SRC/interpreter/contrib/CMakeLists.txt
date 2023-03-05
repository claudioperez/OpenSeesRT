add_library(OpenSeesCommunity OBJECT)
#target_link_libraries(OpenSeesRT PRIVATE OpenSeesCommunity)


add_subdirectory(packages/)
add_subdirectory(uniaxial/)
# add_subdirectory(sys_of_eqn)
# add_subdirectory(packages/YieldSurface/)

target_sources(OpenSeesRT PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/integrator/Newmark/G3Parse_Newmark.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/element/UP-ucsd/TclFourNodeQuadUPCommand.cpp"
    #"element/UWelements/Tcl_generateInterfacePoints.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/adapter/TclActuatorCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/adapter/TclActuatorCorotCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/adapter/TclAdapterCommand.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/element/brick/TclBrickCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/brick/TclTwenty_Node_BrickCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/elasticBeamColumn/TclElasticBeamCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/elastomericBearing/TclElastomericBearingBoucWenCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/elastomericBearing/TclElastomericBearingPlasticityCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/elastomericBearing/TclElastomericBearingUFRPCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/fourNodeQuad/TclFourNodeQuadCommand.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/element/frictionBearing/TclFlatSliderCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/frictionBearing/TclRJWatsonEQSCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/frictionBearing/TclSingleFPCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/frictionBearing/frictionModel/TclModelBuilderFrictionModelCommand.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/element/generic/TclGenericClientCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/generic/TclGenericCopyCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/joint/TclBeamColumnJointCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/joint/TclJoint2dCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/joint/TclJoint3dCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/twoNodeLink/TclTwoNodeLinkCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/zeroLength/TclZeroLength.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/element/nonlinearBeamColumn/updatedLagrangianBeamColumn/TclCyclicModelCommands.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/nonlinearBeamColumn/updatedLagrangianBeamColumn/TclElement2dGNL.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/nonlinearBeamColumn/updatedLagrangianBeamColumn/TclElement2dYS.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/element/nonlinearBeamColumn/forceBeamColumn/TclForceBeamColumnCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/nonlinearBeamColumn/dispBeamColumnInt/TclDispBeamColumnIntCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/nonlinearBeamColumn/beamWithHinges.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/nonlinearBeamColumn/gradientInelasticBeamColumn.cpp"
)


