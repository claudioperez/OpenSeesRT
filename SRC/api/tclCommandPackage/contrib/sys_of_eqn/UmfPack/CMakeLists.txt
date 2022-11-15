target_sources(OpenSeesRT_LinearSystem PRIVATE
  umfpack.cpp
)

find_package(SuiteSparse)
target_link_libraries(OpenSeesRT_LinearSystem PRIVATE -lumfpack)

