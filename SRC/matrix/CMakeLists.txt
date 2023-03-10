#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(matrix OBJECT)

target_sources(OPS_Matrix
    PRIVATE
      eig3.cpp
      Matrix.cpp
      Vector.cpp
      R3vectors.cpp
      ID.cpp
    PUBLIC
      eig3.h
      Matrix.h
      Vector.h
      R3vectors.h
      ID.h
)


target_include_directories(OPS_Matrix PUBLIC ${CMAKE_CURRENT_LIST_DIR})
