#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Matrix
  PRIVATE
    cbdi.cpp
  PUBLIC
    cbdi.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OPS_Matrix PUBLIC ${CMAKE_CURRENT_LIST_DIR})


