#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Parallel_SP
  PRIVATE
    DomainPartitioner.cpp
  PUBLIC
    DomainPartitioner.h
)

target_include_directories(OPS_Parallel_SP PUBLIC ${CMAKE_CURRENT_LIST_DIR})
