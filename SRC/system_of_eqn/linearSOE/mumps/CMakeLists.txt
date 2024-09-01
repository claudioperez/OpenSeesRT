#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Runtime PRIVATE dispatch.cpp)

target_sources(OPS_SysOfEqn
  PRIVATE
#   MumpsParallelSOE.cpp
#   MumpsParallelSolver.cpp
    MumpsSOE.cpp
    MumpsSolver.cpp
  PUBLIC
    MumpsSOE.h
    MumpsSolver.h

)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

