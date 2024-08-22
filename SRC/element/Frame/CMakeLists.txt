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
add_subdirectory(Basic)
add_subdirectory(Local)
add_subdirectory(Total)
add_subdirectory(Elastic)
add_subdirectory(Other/Displ)
add_subdirectory(Other/Force)
add_subdirectory(Other/Mixed)
add_subdirectory(Other/ShearWall)
add_subdirectory(Other/GradientInelastic)
add_subdirectory(Other/UpdatedLagrange)
