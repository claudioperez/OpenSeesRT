#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_include_directories(OPS_Parallel_SP PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(OPS_Parallel_SP
  PRIVATE
    PartitionedDomainEleIter.cpp 
    PartitionedDomainSubIter.cpp
    PartitionedDomain.cpp
  PUBLIC
    PartitionedDomainEleIter.h
    PartitionedDomainSubIter.h
    PartitionedDomain.h
)

