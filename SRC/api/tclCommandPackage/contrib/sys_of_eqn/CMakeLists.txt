# add_subdirectory(UmfPack)
# add_subdirectory(Mumps)
add_subdirectory(PetSc)
add_subdirectory(ITPACK)
add_subdirectory(ThreadedSuperLU)

target_link_libraries(OpenSeesRT PUBLIC OpenSeesRT_LinearSystem)

