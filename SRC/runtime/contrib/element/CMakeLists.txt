
target_sources(OPS_Runtime PRIVATE
  PRIVATE

    "${CMAKE_CURRENT_LIST_DIR}/Frame/beamWithHinges.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/Frame/updatedLagrangianBeamColumn/TclElement2dGNL.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/Frame/updatedLagrangianBeamColumn/TclElement2dYS.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/Frame/updatedLagrangianBeamColumn/TclCyclicModelCommands.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/Other/TclGenericClientCommand.cpp"
)
