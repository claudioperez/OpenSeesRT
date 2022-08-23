# add_library(OpenSeesPyRT SHARED)
# add_library(OpenSeesPyRT OBJECT)

find_package(PythonExtensions REQUIRED)
execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pybind11 --cmakedir OUTPUT_VARIABLE pybind11_ROOT)
find_package(pybind11 CONFIG)

pybind11_add_module(OpenSeesPyRT)
#python_add_library(OpenSeesPyRT)

target_include_directories(OpenSeesPyRT 
  PRIVATE ${Python_INCLUDE_DIRS}
          ${pybind11_INCLUDE_DIRS}
)

set_property(TARGET OpenSeesPyRT PROPERTY POSITION_INDEPENDENT_CODE 1)
target_sources(OpenSeesPyRT PRIVATE
  "OpenSeesPyRT.cpp"
)
# set_target_properties(OpenSeesPyRT PROPERTIES CXX_VISIBILITY_PRESET hidden)
target_link_libraries(OpenSeesPyRT  PRIVATE ${pybind11_LIBRARIES} ${TCL_STUB_LIBRARY})
target_link_libraries(OpenSeesPyRT  PRIVATE OpenSeesRT)

