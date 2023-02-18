
set(PETSC_DIR $ENV{PETSC_DIR})
# set(PETSC_LIBRARIES $ENV{PETSC_DIR}/lib/libpetsc.so)
set(PETSC_LIBRARIES $ENV{CONDA_PREFIX}/lib/libpetsc.so)

message(">>> PETSC: ${PETSC_DIR}")
message(">>> PETSC: ${PETSC_LIBRARIES}")

add_library(PETSc_SOE)

target_link_libraries(PETSc_SOE PRIVATE OPS_SysOfEqn ${PETSC_LIBRARIES})
target_include_directories(PETSc_SOE PRIVATE ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(OpenSeesRT_LinearSystem PRIVATE PETSc_SOE)
target_compile_definitions(OpenSeesRT PRIVATE OPS_PETSC)

target_sources(PETSc_SOE
    PRIVATE
      RT_PetSc.cpp
    #TclPetsc.cpp
#     ActorPetscSOE.cpp
#     badPetscSOE.cpp
      PetscSOE.cpp
      PetscSolver.cpp
      PetscSparseSeqSolver.cpp
#     ShadowPetscSOE.cpp

    PUBLIC
#     ActorPetscSOE.h
      oldpetsc.h
      PetscSOE.h
      PetscSolver.h
      PetscSparseSeqSolver.h
#     ShadowPetscSOE.h
)


#     main.cpp
