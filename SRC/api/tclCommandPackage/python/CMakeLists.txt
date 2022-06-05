add_library(OpenSeesPyRT SHARED)
# add_library(OpenSeesPyRT OBJECT)

find_package(Python)
find_package(pybind11 CONFIG)
target_include_directories(OpenSeesPyRT 
  PRIVATE ${Python_INCLUDE_DIRS}
          ${pybind11_INCLUDE_DIRS}
)
target_sources(OpenSeesPyRT PRIVATE
  "OpenSeesPyRT.cpp"
)
set_target_properties(OpenSeesPyRT PROPERTIES CXX_VISIBILITY_PRESET hidden)
target_link_libraries(OpenSeesPyRT PRIVATE ${pybind11_LIBRARIES} ${TCL_LIBRARIES})

