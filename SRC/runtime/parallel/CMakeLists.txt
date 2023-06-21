
if(MKL_FOUND)
  message(STATUS "MKL was found.")
  set(SCALAPACK_LIBRARIES ${MKL_LIBRARIES})
else()
  message(STATUS "MKL NOT found.")
  find_package(SCALAPACK)
  message(">>> ScaLAPACK: ${ScaLAPACK_LIBRARIES}")
endif()
  
target_link_libraries(OpenSeesRT_Parallel PRIVATE OPS_Parallel OPS_Actor
    ${MPI_CXX_LIBRARIES}
    /usr/lib/libscalapack.so
    ${MPI_Fortran_LIBRARIES}
    ${MPI_CXX_LINK_FLAGS}
)

target_sources(OpenSeesMP PRIVATE mpiMain.cpp)
target_link_libraries(OpenSeesMP PRIVATE OPS_Parallel OPS_Actor  OpenSeesRT ${MPI_LIBRARIES})

target_sources(OpenSeesSP PRIVATE mpiMain.cpp)
target_link_libraries(OpenSeesSP PRIVATE OPS_Parallel OPS_Actor  OpenSeesRT ${MPI_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    /usr/lib/libscalapack.so
    ${MPI_Fortran_LIBRARIES}       
    ${MPI_CXX_LINK_FLAGS}
)
