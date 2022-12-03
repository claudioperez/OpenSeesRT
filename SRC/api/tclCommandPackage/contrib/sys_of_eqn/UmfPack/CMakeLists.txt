
# find_package(SuiteSparse)
target_include_directories(OpenSeesRT_LinearSystem
  PRIVATE "${OPS_SRC_DIR}/system_of_eqn/linearSOE/umfGEN"
)
target_link_libraries(OpenSeesRT_LinearSystem PUBLIC 
  OPS_SysOfEqn_UMF
  /home/claudio/packages/opensees-pypi/src/libg3/OTHER/bin/UMFPACK/libUMFPACK.a
)
target_link_libraries(OpenSeesRT PUBLIC /home/claudio/packages/opensees-pypi/src/libg3/OTHER/bin/UMFPACK/libUMFPACK.a)

target_sources(OpenSeesRT_LinearSystem PRIVATE
  umfpack.cpp
)
