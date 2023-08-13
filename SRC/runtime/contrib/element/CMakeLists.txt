
target_sources(OPS_Runtime PRIVATE
  PRIVATE
#   "Plane/QuadCommands.cpp"
#   "Plane/TclFourNodeQuadUPCommand.cpp"

    "brick/TclBrickCommand.cpp"
    "brick/TclTwenty_Node_BrickCommand.cpp"

    "Frame/TclForceBeamColumnCommand.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/Frame/beamWithHinges.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/Frame/TclDispBeamColumnIntCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/Frame/TclElasticBeamCommand.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/Frame/updatedLagrangianBeamColumn/TclElement2dGNL.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/Frame/updatedLagrangianBeamColumn/TclElement2dYS.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/Frame/updatedLagrangianBeamColumn/TclCyclicModelCommands.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/Other/TclGenericClientCommand.cpp"
)
