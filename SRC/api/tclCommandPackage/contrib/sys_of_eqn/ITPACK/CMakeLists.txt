#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OpenSeesRT_LinearSystem
    PRIVATE
      TclDispatchITPACK.cpp
      ItpackLinSOE.cpp
      ItpackLinSolver.cpp
      ItpackLinSOE.h
      ItpackLinSolver.h
    PUBLIC

)

add_subdirectory(ITPACK)
target_link_libraries(OpenSeesRT_LinearSystem PRIVATE ITPACK)
target_include_directories(OpenSeesRT_LinearSystem PUBLIC ${CMAKE_CURRENT_LIST_DIR})
