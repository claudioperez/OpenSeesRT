
add_library(SuperLU_SOE)

target_sources(SuperLU_SOE
  PRIVATE
    ThreadedSuperLU.cpp
# PUBLIC
#   ThreadedSuperLU.h
)


add_subdirectory(SuperLU_MT-3/SRC)

target_link_libraries(SuperLU_SOE PRIVATE OPS_SysOfEqn SuperLU_MT)

target_link_libraries(OpenSeesRT  PRIVATE SuperLU_SOE)

