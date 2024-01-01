#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#
# Claudio M. Perez
#
add_library(OPS_Quadrature)

target_sources(OPS_Quadrature
  PRIVATE
    GaussNURBS.cpp 
  PUBLIC
    GaussNURBS.h
)

target_link_libraries(OPS_Element PRIVATE OPS_Quadrature)

add_subdirectory(Frame)

