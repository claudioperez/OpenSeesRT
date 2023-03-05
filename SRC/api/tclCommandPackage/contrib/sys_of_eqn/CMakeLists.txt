
add_library(OpenSeesRT_LinearSystem OBJECT)
target_compile_definitions(OpenSeesRT_LinearSystem PUBLIC USE_TCL_STUBS)
target_link_libraries(OpenSeesRT_LinearSystem PUBLIC ${TCL_LIBRARY})
target_link_libraries(OpenSeesRT PUBLIC OpenSeesRT_LinearSystem)

# add_subdirectory(UmfPack)
# add_subdirectory(Mumps)
add_subdirectory(PetSc)
add_subdirectory(ITPACK)
add_subdirectory(ThreadedSuperLU)

