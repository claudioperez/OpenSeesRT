#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(coordTransformation "")

target_link_libraries(coordTransformation PRIVATE OPS_Recorder)

target_sources(coordTransformation
  PRIVATE
    CrdTransf.cpp
    LinearCrdTransf2d.cpp
    LinearCrdTransf2d02.cpp
    LinearCrdTransf3d.cpp
    LinearCrdTransf3d02.cpp
    PDeltaCrdTransf2d.cpp
    PDeltaCrdTransf3d.cpp
    CorotCrdTransf2d.cpp
    CorotCrdTransf3d.cpp
    CorotCrdTransf3d02.cpp
    CorotCrdTransfWarping2d.cpp
  PUBLIC
    CrdTransf.h
    LinearCrdTransf2d.h
    LinearCrdTransf2d02.h
    LinearCrdTransf3d.h
    LinearCrdTransf3d02.h
    PDeltaCrdTransf2d.h
    PDeltaCrdTransf3d.h
    CorotCrdTransf2d.h
    CorotCrdTransf3d.h
    CorotCrdTransf3d02.h
    CorotCrdTransfWarping2d.h

    ASDEICR.h
    ASDShellQ4CorotationalTransformation.h
    ASDShellQ4LocalCoordinateSystem.h
    ASDShellQ4Transformation.h
)

target_include_directories(coordTransformation PUBLIC ${CMAKE_CURRENT_LIST_DIR})
