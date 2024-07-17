#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Transform
    PRIVATE
        CorotFrameTransf3d.cpp
        LinearFrameTransf3d.cpp
        PDeltaFrameTransf3d.cpp
    PUBLIC
        CorotFrameTransf3d.h
        LinearFrameTransf3d.h
        PDeltaFrameTransf3d.h
)
target_include_directories(OPS_Transform PUBLIC ${CMAKE_CURRENT_LIST_DIR})

