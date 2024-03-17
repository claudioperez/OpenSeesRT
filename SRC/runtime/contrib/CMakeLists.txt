add_library(OpenSeesCommunity OBJECT)



add_subdirectory(packages/)
add_subdirectory(uniaxial/)
add_subdirectory(element/)
# add_subdirectory(packages/YieldSurface/)

target_sources(OPS_Runtime PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/integrator/Newmark/G3Parse_Newmark.cpp"

    #"element/UWelements/Tcl_generateInterfacePoints.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/adapter/TclActuatorCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/adapter/TclActuatorCorotCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/adapter/TclAdapterCommand.cpp"


    "${CMAKE_CURRENT_LIST_DIR}/element/frictionBearing/TclFlatSliderCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/frictionBearing/TclRJWatsonEQSCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/frictionBearing/TclSingleFPCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/element/frictionBearing/frictionModel/TclModelBuilderFrictionModelCommand.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/element/zeroLength/TclZeroLength.cpp"
)


