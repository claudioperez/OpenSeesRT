#------------------------------------------------------------------------------
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#------------------------------------------------------------------------------

target_sources(OPS_Element
  PRIVATE
    PrismFrame2d.cpp
    PrismFrame3d.cpp

    CubicFrame3d.cpp

    MixedFrame3d.cpp

    ForceFrame3d.cpp
    ForceDeltaFrame3d.cpp

    EulerFrame3d.cpp
    EulerDeltaFrame3d.cpp

  PUBLIC
    PrismFrame2d.h
    PrismFrame3d.h
    CubicFrame3d.h
    ExactFrame3d.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

