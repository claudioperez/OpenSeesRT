#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

add_executable(OpenSees
        tclMain.cpp 
#       tclAppInit.cpp
)

target_link_libraries(OpenSees PRIVATE ${TCL_LIBRARY} OpenSeesRT)

