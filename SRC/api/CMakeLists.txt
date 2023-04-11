#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Api OBJECT)

target_sources(OPS_Api
    PUBLIC
    # packages.cpp
    PUBLIC
    # elementAPI.h
    # packages.h
)
# target_include_directories(OPS_Api PUBLIC ${CMAKE_CURRENT_LIST_DIR})

# if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/tclCommandPackage/CMakeLists.txt")
# add_subdirectory(tclCommandPackage)
# endif()

