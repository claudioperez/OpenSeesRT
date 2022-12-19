# cmake_minimum_required(VERSION 3.24)
# project(UCSD)
add_subdirectory(nDarray)

add_library(UCSD_FiniteDeformation)
target_include_directories(UCSD_FiniteDeformation PRIVATE "${CMAKE_CURRENT_LIST_DIR}/nDarray/")
add_subdirectory(FiniteDeformation)

