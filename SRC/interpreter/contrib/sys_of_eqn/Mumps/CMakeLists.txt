

target_compile_definitions(OpenSeesRT_LinearSystem PRIVATE _MUMPS)
# find_package(BLAS)
# find_package(LAPACK)
set(SCALAPACK_LIBRARY "/home/claudio/mambaforge/envs/py39/lib/libscalapack.so")
# set(SCALAPACK_links true)
# set(SCALAPACK_INCLUDE_DIR ""
# find_package(SCALAPACK)

# --- Fetch and build MUMPS
include(FetchContent)
set(CMAKE_TLS_VERIFY true)
set(FETCHCONTENT_UPDATES_DISCONNECTED_MUMPSUPSTREAM true)

find_package(MPI COMPONENTS C Fortran REQUIRED)

# FetchContent_Declare(ScaLAPACKupstream
#  GIT_REPOSITORY https://github.com/scivision/scalapack.git
#  GIT_TAG v2.2.1.4
# )

# FetchContent_MakeAvailable(ScaLAPACKupstream)
# list( APPEND CMAKE_PREFIX_PATH ${FETCHCONTENT_BASE_DIR}/scalapackupstream-subbuild )

# set (SCALAPACK_link ScaLAPACK)

# find_package( SCALAPACK REQUIRED HINTS "${SCALAPACK_SOURCE_DIR}" "${ScaLAPACK_DIR}")

# add_subdirectory("${FETCHCONTENT_BASE_DIR}/scalapackupstream-subbuild/CMakeFiles/")
if( FALSE AND NOT TARGET ScaLAPACK )
  # Externally installed scalapack
  find_package( ScaLAPACK QUIET )

  if( TARGET ScaLAPACK::ScaLAPACK )
    message( STATUS "Found external ScaLAPACK installation" )
    add_library( ChronusQ::ScaLAPACK ALIAS ScaLAPACK::ScaLAPACK )
  else()
    list( APPEND CMAKE_PREFIX_PATH ${FETCHCONTENT_BASE_DIR}/scalapack-build )
    find_package( ScaLAPACK QUIET )

    if( TARGET ScaLAPACK::ScaLAPACK )
      message( STATUS "Found previously installed ScaLAPACK" )
      add_library( ChronusQ::ScaLAPACK ALIAS ScaLAPACK::ScaLAPACK )

    else()
      message( STATUS "Building ScaLAPACK" )
      FetchContent_Declare( ReferenceScaLAPACK
        GIT_REPOSITORY https://github.com/Reference-ScaLAPACK/scalapack.git
        GIT_TAG bc6cad585362aa58e05186bb85d4b619080c45a9
      )
      FetchContent_MakeAvailable( ReferenceScaLAPACK )
      set_property( TARGET scalapack APPEND PROPERTY COMPILE_FLAGS "-Wno-implicit-function-declaration" )
      set_property( TARGET scalapack APPEND PROPERTY INTERFACE_COMPILE_FLAGS "-Wno-implicit-function-declaration" )
      add_library( ScaLAPACK ALIAS scalapack )
    endif()
  endif()
endif()

# add_library( ScaLAPACK ALIAS scalapack )
# find_package( scalapack )
message(">>> SCALAPACK:         ${ScaLAPACK}")
message(">>> SCALAPACK_LIBRARY: ${SCALAPACK_LIBRARY}")
message(">>> SCALAPACK_LIBRARY: ${SCALAPACK_LIBRARIES}")
message(">>> SCALAPACK_links:   ${SCALAPACK_links}")
# message(">>> ScaLAPACK:         ${ScaLAPACK::ScaLAPACK}")
message(">>> ScaLAPACKupstream: ${ScaLAPACK_LIBRARIES}")
message(">>> : ${scalapack_LIBRARIES}")
message(">>> scalapack:         ${scalapack}")
set (SCALAPACK_LIBRARY scalapack CACHE FILEPATH "scalapack")
set (SCALAPACK_link scalapack CACHE FILEPATH "scalapack")
#message(">>> SCALAPACK:         ${SCALAPACK::SCALAPACK}")


FetchContent_Declare(MUMPSupstream
  GIT_REPOSITORY https://github.com/scivision/mumps.git
  GIT_TAG v5.5.1.7
)

FetchContent_MakeAvailable(MUMPSupstream)
#-------------------------------------------------------------


target_sources(OpenSeesRT_LinearSystem PRIVATE
  TclDispatchMumps.cpp
  "${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSOE.cpp"
  "${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSolver.cpp"
)

target_include_directories(OpenSeesRT_LinearSystem PRIVATE
  "${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/"
)

# target_link_libraries( PRIVATE -lmumps)
target_link_libraries(OpenSeesRT_LinearSystem PRIVATE MUMPS::MUMPS MPI::MPI_C)


