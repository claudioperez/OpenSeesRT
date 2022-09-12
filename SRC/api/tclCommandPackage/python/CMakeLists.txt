# add_library(OpenSeesPyRT SHARED)
# add_library(OpenSeesPyRT OBJECT)

find_package(PythonExtensions REQUIRED)
find_package(pybind11 CONFIG)

pybind11_add_module(OpenSeesPyRT)

target_include_directories(OpenSeesPyRT 
  PRIVATE ${Python_INCLUDE_DIRS}
          ${pybind11_INCLUDE_DIRS}
)

target_sources(OpenSeesPyRT PRIVATE
  "OpenSeesPyRT.cpp"
)
# set_target_properties(OpenSeesPyRT PROPERTIES CXX_VISIBILITY_PRESET hidden)
set_property(TARGET OpenSeesPyRT PROPERTY POSITION_INDEPENDENT_CODE 1)
target_link_libraries(OpenSeesPyRT  PRIVATE ${pybind11_LIBRARIES} ${TCL_STUB_LIBRARY})
target_link_libraries(OpenSeesPyRT  PRIVATE OpenSeesRT)

