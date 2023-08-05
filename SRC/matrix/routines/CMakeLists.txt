#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Matrix
  PRIVATE
    eig3.cpp
  PUBLIC
    eig3.h
)

target_include_directories(OPS_Matrix PUBLIC ${CMAKE_CURRENT_LIST_DIR})
