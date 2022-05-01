add_library(OpenSeesPyRT SHARED)

find_package(Python)
find_package(pybind11 CONFIG)
target_include_directories(OpenSeesPyRT 
  PRIVATE ${Python_INCLUDE_DIRS}
          ${pybind11_INCLUDE_DIRS}
)
target_sources(OpenSeesPyRT PRIVATE
  "OpenSeesPyRT.cpp"
)
target_link_libraries(OpenSeesPyRT PRIVATE ${pybind11_LIBRARIES} ${TCL_LIBRARIES})

