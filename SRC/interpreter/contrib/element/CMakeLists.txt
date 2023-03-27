
target_sources(OpenSeesRT PRIVATE
  PRIVATE
    "Plane/QuadCommands.cpp"
    "Plane/TclFourNodeQuadUPCommand.cpp"

    "brick/TclBrickCommand.cpp"
    "brick/TclTwenty_Node_BrickCommand.cpp"

    "Frame/TclForceBeamColumnCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/Frame/updatedLagrangianBeamColumn/TclElement2dGNL.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/Frame/updatedLagrangianBeamColumn/TclElement2dYS.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/Frame/updatedLagrangianBeamColumn/TclCyclicModelCommands.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/Frame/TclElasticBeamCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/Frame/TclDispBeamColumnIntCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/Frame/beamWithHinges.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/Frame/gradientInelasticBeamColumn.cpp"

    "Shell/IO_ASDShellQ4.cpp"
    "Shell/IO_ShellANDeS.cpp"
    "Shell/IO_ShellDKGQ.cpp"
    "Shell/IO_ShellDKGT.cpp"
    "Shell/IO_ShellMITC4.cpp"
    "Shell/IO_ShellMITC4Thermal.cpp"
    "Shell/IO_ShellMITC9.cpp"
    "Shell/IO_ShellNLDKGQ.cpp"
    "Shell/IO_ShellNLDKGQThermal.cpp"
    "Shell/IO_ShellNLDKGT.cpp"


    "${CMAKE_CURRENT_LIST_DIR}/Other/TclGenericClientCommand.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/Other/TclGenericCopyCommand.cpp"
)
