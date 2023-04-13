#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Parallel
  PRIVATE
    DomainPartitioner.cpp
  PUBLIC
    DomainPartitioner.h
)

target_include_directories(OPS_Parallel PUBLIC ${CMAKE_CURRENT_LIST_DIR})
