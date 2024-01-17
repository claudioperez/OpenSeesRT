# Building OpenSeesRT on TACC
#
# Modify `CMakeLists.txt` files:
# 
#   `/src/libg3/CMakeLists.txt`:
#   Uncomment `find_package(MKL)`
#   `/src/libg3/SRC/runtime/parallel/CMakeLists.txt`:
#   Remove the two `/usr/lib/libscalapack.so` lines
#
#   Also modify `setup.py` to specify `"-j56"` towards the end of the file,
#   utilizing all cores for parallel compilation.
# 
#  Links:
#  [1] Good Conduct on Frontera: <https://frontera-portal.tacc.utexas.edu/user-guide/conduct/>
#  [2] Using Modules to Manage your Environment: 
#      <https://frontera-portal.tacc.utexas.edu/user-guide/admin/#using-modules-to-manage-your-environment>
#  [3] Architecture-Specific Flags: 
#      <https://frontera-portal.tacc.utexas.edu/user-guide/building/#architecture-specific-flags>


#
#.  Load the appropriate modules. See [2].
#
    module load impi
    module load python3/3.9.2  # this is currently the latest available python 3
#
#.  Configure environment variables
#
    export CC=$(which mpicc)
    export CXX=$(which mpicxx)
    export FC=$(which mpif90)
    export LD_PRELOAD=/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64_lin/libmkl_def.so:/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64_lin/libmkl_avx2.so:/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64_lin/libmkl_core.so:/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64_lin/libmkl_intel_lp64.so:/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64_lin/libmkl_intel_thread.so:/opt/intel/compilers_and_libraries_2020.1.217/linux/compiler/lib/intel64_lin/libiomp5.so

#
#.  Install `OpenSeesRT`
#
    
    # preserving the build directory
#   python3 setup.py cmake
#   make -C build/temp.linux-x86_64-cpython-39_local/ OpenSeesRT -j56
#   export OPENSEESRT_LIB=/work2/07506/usr83847/frontera/OpenSeesRT/build/temp.linux-x86_64-cpython-39_local/src/libg3/SRC/runtime/libOpenSeesRT.so
    
    # or not
#   pip3 install --user -e .  # installs to $HOME/.local

#
#.  Run
#
# Before running any code that requires `OpenSeesRT`, the environment setup
# should be configured as follows:

#   module load impi
#   module load python3/3.9.2
#   export LD_PRELOAD=/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64_lin/libmkl_def.so:/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64_lin/libmkl_avx2.so:/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64_lin/libmkl_core.so:/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64_lin/libmkl_intel_lp64.so:/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64_lin/libmkl_intel_thread.so:/opt/intel/compilers_and_libraries_2020.1.217/linux/compiler/lib/intel64_lin/libiomp5.so
#   
#   export PYTHONPATH=$PYTHONPATH:$(pwd)

# Optional: Recompile specifying the instruction set. See [3].
# 
# Edit `src/libg3/CMakeLists.txt` and append the following optimization options:
# 
#     $<$<CONFIG:RELEASE>:-mkl=cluster>
#     $<$<CONFIG:RELEASE>:-xCORE-AVX512>
# 
#     python3 setup.py cmake
#     cd build/temp.linux-x86_64-cpython-39_pypa
#     make OpenSeesRT -j52
#     
#     export OPENSEESRT_LIB=$SCRATCH/OpenSeesRT/build/temp.linux-x86_64-cpython-39_pypa/src/libg3/SRC/runtime/libOpenSeesRT.so

