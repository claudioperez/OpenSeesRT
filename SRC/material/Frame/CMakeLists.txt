#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
  PRIVATE
    ElasticIsotropicBeamFiber2d.cpp
    ElasticIsotropicBeamFiber.cpp
    J2BeamFiber2d.cpp
    J2BeamFiber3d.cpp
    BeamFiberMaterial.cpp
    BeamFiberMaterial2d.cpp
    BeamFiberMaterial2dPS.cpp
  PUBLIC
    J2BeamFiber2d.h
    J2BeamFiber3d.h
    BeamFiberMaterial.h
    BeamFiberMaterial2d.h
    BeamFiberMaterial2dPS.h
    ElasticIsotropicBeamFiber2d.h
    ElasticIsotropicBeamFiber.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
