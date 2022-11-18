
# find_package(SuiteSparse)
target_include_directories(OpenSeesRT_LinearSystem
  PRIVATE "${OPS_SRC_DIR}/system_of_eqn/linearSOE/umfGEN"
)
target_link_libraries(OpenSeesRT_LinearSystem PRIVATE -lumfpack)

target_sources(OpenSeesRT_LinearSystem PRIVATE
  umfpack.cpp
)
