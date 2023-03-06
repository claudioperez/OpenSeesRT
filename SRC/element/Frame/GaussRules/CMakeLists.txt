#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_include_directories(OPS_Element PUBLIC $(CMAKE_CURRENT_LIST_DIR))

target_sources(OPS_Element
  PRIVATE
        BeamIntegration.cpp
        LobattoBeamIntegration.cpp
        LegendreBeamIntegration.cpp
        RadauBeamIntegration.cpp
        NewtonCotesBeamIntegration.cpp
        HingeMidpointBeamIntegration.cpp
        HingeEndpointBeamIntegration.cpp
        HingeRadauBeamIntegration.cpp
        HingeRadauTwoBeamIntegration.cpp
        UserDefinedBeamIntegration.cpp
        SimpsonBeamIntegration.cpp
        UserDefinedHingeIntegration.cpp
        DistHingeIntegration.cpp
        RegularizedHingeIntegration.cpp
        TrapezoidalBeamIntegration.cpp
        CompositeSimpsonBeamIntegration.cpp
        FixedLocationBeamIntegration.cpp
        LowOrderBeamIntegration.cpp
        MidDistanceBeamIntegration.cpp
        #GaussQBeamIntegration.cpp
	ChebyshevBeamIntegration.cpp
  PUBLIC
        BeamIntegration.h
        LobattoBeamIntegration.h
        LegendreBeamIntegration.h
        RadauBeamIntegration.h
        NewtonCotesBeamIntegration.h
        HingeMidpointBeamIntegration.h
        HingeEndpointBeamIntegration.h
        HingeRadauBeamIntegration.h
        HingeRadauTwoBeamIntegration.h
        UserDefinedBeamIntegration.h
        SimpsonBeamIntegration.h
        UserDefinedHingeIntegration.h
        DistHingeIntegration.h
        RegularizedHingeIntegration.h
        TrapezoidalBeamIntegration.h
        CompositeSimpsonBeamIntegration.h
        FixedLocationBeamIntegration.h
        LowOrderBeamIntegration.h
        MidDistanceBeamIntegration.h
        #GaussQBeamIntegration.h
	ChebyshevBeamIntegration.h
)
