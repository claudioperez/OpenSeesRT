target_sources(OPS_Runtime PRIVATE commands.cpp)

target_sources(OPS_Domain
  PRIVATE
    PeriDomain.cpp
  PUBLIC
    PeriDomain.h
)
target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})
