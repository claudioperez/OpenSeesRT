target_sources(OPS_Domain
  PRIVATE
    ParticleDomain.cpp
  PUBLIC
    ParticleDomain.h
)
target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})
