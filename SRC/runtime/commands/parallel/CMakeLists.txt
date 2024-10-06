#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
find_package(MPI)
if (NOT MPI_FOUND)
  return()
endif()

add_library(OpenSeesRT_Parallel OBJECT)

target_include_directories(OpenSeesRT_Parallel
  PRIVATE
  ${OPS_SRC_DIR}/domain/partitioner/
  ${OPS_SRC_DIR}/analysis/analysis/
)

target_compile_definitions(OpenSeesRT_Parallel PRIVATE 
  _PARALLEL_INTERPRETERS
  _PARALLEL_PROCESSING     # OpenSeesSP
)
target_link_libraries(OpenSeesRT_Parallel PRIVATE 
    OPS_Runtime 
    OPS_Parallel 
    OPS_Actor ${MPI_CXX_LIBRARIES})

target_sources(OpenSeesRT_Parallel PRIVATE
  machine.cpp
)

add_library(LibOpenSeesSP  SHARED  EXCLUDE_FROM_ALL )
set_property(TARGET LibOpenSeesSP PROPERTY POSITION_INDEPENDENT_CODE 1)

add_library(LibOpenSeesMP  SHARED  EXCLUDE_FROM_ALL )
set_property(TARGET LibOpenSeesMP PROPERTY POSITION_INDEPENDENT_CODE 1)

if(MKL_FOUND)
  message(STATUS "MKL was found.")
  set(SCALAPACK_LIBRARIES ${MKL_LIBRARIES})
else()
  message(STATUS "MKL NOT found.")
  find_package(SCALAPACK)
  message(">>> ScaLAPACK: ${ScaLAPACK_LIBRARIES}")
endif()
  
target_link_libraries(OpenSeesRT_Parallel PRIVATE 
    OPS_Parallel OPS_Actor
    ${MPI_CXX_LIBRARIES}
#   /usr/lib/libscalapack.so
    ${MPI_Fortran_LIBRARIES}
    ${MPI_CXX_LINK_FLAGS}
)

target_link_libraries(LibOpenSeesMP PRIVATE 
        OpenSeesRT_Parallel
	OPS_Parallel 
	OPS_Actor  
	OpenSeesRT 
	OPS_Runtime 
	${MPI_LIBRARIES})

target_sources(LibOpenSeesMP PRIVATE 
    communicate.cpp
    ${OPS_SRC_DIR}/parallel/OpenSeesMP.cpp
)

target_sources(LibOpenSeesSP PRIVATE 
    partition.cpp
    ${OPS_SRC_DIR}/parallel/OpenSeesSP.cpp
)

target_link_libraries(LibOpenSeesSP 
  PRIVATE 
    OpenSeesRT_Parallel
    OPS_Parallel 
    OPS_Partition
    OPS_Analysis
    OPS_Runtime 
    OPS_Actor  OpenSeesRT ${MPI_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    METIS
#   /usr/lib/libscalapack.so
    ${MPI_Fortran_LIBRARIES}       
    ${MPI_CXX_LINK_FLAGS}
)
