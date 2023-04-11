add_library(OpenSeesCommunity OBJECT)
#target_link_libraries(OpenSeesRT PRIVATE OpenSeesCommunity)


add_subdirectory(packages/)
add_subdirectory(uniaxial/)
add_subdirectory(element/)
add_subdirectory(sys_of_eqn)
# add_subdirectory(packages/YieldSurface/)

target_sources(OpenSeesRT PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/integrator/Newmark/G3Parse_Newmark.cpp"

    #"element/UWelements/Tcl_generateInterfacePoints.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/adapter/TclActuatorCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/adapter/TclActuatorCorotCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/adapter/TclAdapterCommand.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/element/elastomericBearing/TclElastomericBearingBoucWenCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/elastomericBearing/TclElastomericBearingPlasticityCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/elastomericBearing/TclElastomericBearingUFRPCommand.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/element/frictionBearing/TclFlatSliderCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/frictionBearing/TclRJWatsonEQSCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/frictionBearing/TclSingleFPCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/frictionBearing/frictionModel/TclModelBuilderFrictionModelCommand.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/element/joint/TclBeamColumnJointCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/joint/TclJoint2dCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/joint/TclJoint3dCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/twoNodeLink/TclTwoNodeLinkCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/zeroLength/TclZeroLength.cpp"
)


