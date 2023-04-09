#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
  PUBLIC
    BeamGT.cpp
    MasonPan12.cpp
    MasonPan3D.cpp
  PRIVATE
    BeamGT.h
    MasonPan12.h
    MasonPan3D.h
)    

if (UNIX)
  set_source_files_properties(
      BeamGT.cpp
      MasonPan12.cpp
      MasonPan3D.cpp
      APPEND PROPERTY CMAKE_CXX_FLAGS -w
  )
endif()
