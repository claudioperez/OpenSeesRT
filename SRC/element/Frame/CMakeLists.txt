#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Runtime PRIVATE commands.cpp)

target_sources(OPS_Element
  PRIVATE
    BasicFrame3d.cpp

    PrismFrame2d.cpp
    PrismFrame3d.cpp

    CubicFrame3d.cpp

    MixedFrame3d.cpp

    ForceFrame3d.cpp
    ForceDeltaFrame3d.cpp

    EulerFrame3d.cpp
    EulerDeltaFrame3d.cpp
  PUBLIC
    BasicFrame3d.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

# add_subdirectory(Basic)
add_subdirectory(Elastic)
add_subdirectory(Other/Displ)
add_subdirectory(Other/Force)
add_subdirectory(Other/Mixed)
add_subdirectory(Other/ShearWall)
add_subdirectory(Other/GradientInelastic)
add_subdirectory(Other/UpdatedLagrange)
