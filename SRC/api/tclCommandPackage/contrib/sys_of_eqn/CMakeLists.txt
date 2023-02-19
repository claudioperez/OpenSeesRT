# add_subdirectory(UmfPack)
# add_subdirectory(Mumps)
add_subdirectory(ITPACK)

target_link_libraries(OpenSeesRT PUBLIC OpenSeesRT_LinearSystem)

