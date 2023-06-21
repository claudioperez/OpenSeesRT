#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#                            External Libraries
# - BLAS_LIBRARIES
# - BLAS_INCLUDE_DIRS
#
# - LAPACK_LIBRARIES
# - LAPACK_INCLUDE_DIRS
#
# - ARPACK_LIBRARIES
#
# - SUPERLU_LIBRARIES
# - SUPERLU_INCLUDE_DIRS
#
#
#==============================================================================
# Synopsis
# - opensees_load(<PACKAGE> [FLAGS] [<PATHS>])
#
# Flags:
# - FIND:   Use CMake to find library on system, fail if not found
#
# Keyword arguments:
#   Provide specific paths for library.
#
# - BUNDLED <path/to/OTHER/LIB/>
#   Provide path to OpenSees bundled library directory containing a 
#   CMakeLists.txt file.
#
# - LIBRARY <path/to/lib.a> INCLUDE <path/to/include/>
#
# - CONAN <conan-package/version>
#   Point to a conan package.
#
#----------------------------------------------------------------
# set(TCL_LIBRARIES /home/claudio/miniconda3/envs/elle/lib/libtclstub8.6.a)
# set(TCL_INCLUDE_PATH /home/claudio/miniconda3/envs/elle/include/)
# set(TCL_STUB_LIBRARY ${TCL_LIBRARIES})


opensees_load(BLAS                                           #FIND)
  #LIBRARY /home/claudio/miniconda3/envs/intel/lib/libmkl_rt.so
  LIBRARY /usr/lib/libblas.so.3
)

opensees_load(LAPACK                                         #FIND)
  #LIBRARY /home/claudio/miniconda3/envs/intel/lib/libmkl_rt.so
  LIBRARY /usr/lib/liblapack.so.3
)

opensees_load(METIS                                          FIND)

opensees_load(HDF5                                           FIND)

opensees_load(MySQL                                          FIND)

# find_package(Python COMPONENTS Development)

# Integrated exteral libraries
opensees_load(FEDEAS_Uniaxial LIBRARY FALSE)


