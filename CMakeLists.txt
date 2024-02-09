cmake_minimum_required(VERSION 3.16)
project(OpenSeesRT)

add_subdirectory(src/libg3)

if (TARGET OpenSeesPyRT)
  install(TARGETS 
      OpenSeesRT OpenSeesPyRT LIBRARY DESTINATION .)
else()
  install(TARGETS 
      OpenSeesRT LIBRARY DESTINATION .)
endif()
