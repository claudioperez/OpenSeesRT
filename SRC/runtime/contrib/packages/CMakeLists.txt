#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

add_subdirectory(YieldSurface)

if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/UCSD_CompGeoMech/CMakeLists.txt)
  add_subdirectory(UCSD_CompGeoMech)
endif ()
