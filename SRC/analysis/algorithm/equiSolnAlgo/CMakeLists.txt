#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Analysis
    PRIVATE
      EquiSolnAlgo.cpp 
      ExpressNewton.cpp
      Linear.cpp 
      NewtonRaphson.cpp
      ModifiedNewton.cpp 
      NewtonLineSearch.cpp 
      Broyden.cpp 
      BFGS.cpp
      KrylovNewton.cpp 
      PeriodicNewton.cpp 
      AcceleratedNewton.cpp
      NewtonHallM.cpp
    PUBLIC
      EquiSolnAlgo.h 
      ExpressNewton.h
      Linear.h 
      NewtonRaphson.h
      ModifiedNewton.h 
      NewtonLineSearch.h 
      Broyden.h 
      BFGS.h
      KrylovNewton.h 
      PeriodicNewton.h 
      AcceleratedNewton.h
      NewtonHallM.h
)

target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(accelerator)
add_subdirectory(search)

