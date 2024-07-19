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
  PUBLIC
    BasicFrame3d.h
)
add_subdirectory(Elastic)
add_subdirectory(ForceInterp)
add_subdirectory(Mixed)
add_subdirectory(Basic)
add_subdirectory(Displ)
add_subdirectory(DisplInterpShear)
add_subdirectory(GradientInelastic)
add_subdirectory(UpdatedLagrange)
