#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

add_executable(OpenSees tclMain.cpp)

if (WIN32)
  target_link_libraries(OpenSees PRIVATE Netapi32.lib)
  target_link_libraries(OpenSees PRIVATE ${TCL_LIBRARIES} OpenSeesRT)
else()
  target_link_libraries(OpenSees PRIVATE ${TCL_LIBRARY} OpenSeesRT)
endif()

# if(NOT EXISTS "${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
# 
#   FILE(COPY ${CONAN_LIB_DIRS_TCL}/tcl8.6
#        DESTINATION ${PROJECT_BINARY_DIR}/lib/tcl8.6
#        FILES_MATCHING PATTERN .tcl)
#      
# else ()
#      FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/tcl8.6)
#      add_custom_command(
#       TARGET OpenSees POST_BUILD
#       COMMENT("Copying init.tcl to ${CMAKE_CURRENT_BINARY_DIR}")
#       COMMAND ${CMAKE_COMMAND} -E copy
#         ${CONAN_LIB_DIRS_TCL}/tcl8.6/init.tcl
#         ${CMAKE_CURRENT_BINARY_DIR}/lib/tcl8.6
#      )
# endif()
