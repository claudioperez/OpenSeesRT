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
      ForceFrame3d.cpp
      ForceFrame3d01.cpp
      ForceDeltaFrame3d.cpp
      DisplEulerFrame3d.cpp
    PUBLIC
      PrismFrame3d.h
      PrismFrame2d.h
      CubicFrame3d.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

