#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Analysis
    PRIVATE
#      ExplicitDifference.cpp
#      ArcLength.cpp
#      ArcLength1.cpp
#      BackwardEuler.cpp
#      CentralDifference.cpp
       EigenIntegrator.cpp                  # Jun Peng
       IncrementalIntegrator.cpp
       Integrator.cpp
#      LoadControl.cpp
#      HarmonicSteadyState.cpp
#      LoadPath.cpp
#      MinUnbalDispNorm.cpp
#      Newmark.cpp
#      Newmark1.cpp
       StaticIntegrator.cpp
#      StagedLoadControl.cpp
       TransientIntegrator.cpp
    PUBLIC
#      ExplicitDifference.h
#      ArcLength.h
#      ArcLength1.h
#      BackwardEuler.h
#      CentralDifference.h
       EigenIntegrator.h
       IncrementalIntegrator.h
       Integrator.h
#      LoadControl.h
#      HarmonicSteadyState.h
#      LoadPath.h
#      MinUnbalDispNorm.h
#      Newmark.h
#      Newmark1.h
       StaticIntegrator.h
#      StagedLoadControl.h
       TransientIntegrator.h
)

# target_sources(OPS_Analysis
#     PRIVATE
#        Collocation.cpp                      # Andreas Schellenberg
#        AlphaOS.cpp
#        AlphaOS_TP.cpp
#        AlphaOSGeneralized.cpp
#        AlphaOSGeneralized_TP.cpp            # Andreas Schellenberg
#        CentralDifferenceAlternative.cpp     # fmk
#        CentralDifferenceNoDamping.cpp
#        DisplacementControl.cpp
#        DistributedDisplacementControl.cpp
#        EQPath.cpp                           # M. Salehi
#        GeneralizedAlpha.cpp
#        HHT.cpp
#        HHT_TP.cpp
#        HHTExplicit.cpp
#        HHTExplicit_TP.cpp
#        HHTGeneralized.cpp
#        HHTGeneralized_TP.cpp
#        HHTGeneralizedExplicit.cpp
#        Houbolt.cpp                      # krm
#        HSConstraint.cpp
#        KRAlphaExplicit.cpp
#        KRAlphaExplicit_TP.cpp
#        NewmarkExplicit.cpp              # Andreas Schellenberg
#        ParkLMS3.cpp
#        TRBDF2.cpp
#        TRBDF3.cpp
#        WilsonTheta.cpp
# #      StagedNewmark.cpp
#        GimmeMCK.cpp                     # MHS
#     PUBLIC
#        Collocation.h
#        AlphaOS.h
#        AlphaOS_TP.h
#        AlphaOSGeneralized.h
#        AlphaOSGeneralized_TP.h
#        CentralDifferenceAlternative.h
#        CentralDifferenceNoDamping.h
#        DisplacementControl.h
#        DistributedDisplacementControl.h
#        EQPath.h
#        GeneralizedAlpha.h
#        HHT.h
#        HHT_TP.h
#        HHTExplicit.h
#        HHTExplicit_TP.h
#        HHTGeneralized.h
#        HHTGeneralized_TP.h
#        HHTGeneralizedExplicit.h
#        Houbolt.h
#        HSConstraint.h
#        KRAlphaExplicit.h
#        KRAlphaExplicit_TP.h
#        NewmarkExplicit.h
#        ParkLMS3.h
#        TRBDF2.h
#        TRBDF3.h
#        WilsonTheta.h
# #      StagedNewmark.h
#        GimmeMCK.h
# )
# 
# target_sources(OPS_Analysis
#     PRIVATE
#        CollocationHSFixedNumIter.cpp
#        CollocationHSIncrLimit.cpp
#        CollocationHSIncrReduct.cpp
#        HHTHSFixedNumIter.cpp
#        HHTHSFixedNumIter_TP.cpp
#        HHTHSIncrLimit.cpp
#        HHTHSIncrLimit_TP.cpp
#        HHTHSIncrReduct.cpp
#        HHTHSIncrReduct_TP.cpp
#        HHTGeneralizedExplicit_TP.cpp
#        NewmarkHSFixedNumIter.cpp        # Andreas Schellenberg
#        NewmarkHSIncrLimit.cpp           # Andreas Schellenberg
#        NewmarkHSIncrReduct.cpp          # Andreas Schellenberg
#     PUBLIC
#        CollocationHSFixedNumIter.h
#        CollocationHSIncrLimit.h
#        CollocationHSIncrReduct.h
#        HHTHSFixedNumIter.h
#        HHTHSFixedNumIter_TP.h
#        HHTHSIncrLimit.h
#        HHTHSIncrLimit_TP.h
#        HHTHSIncrReduct.h
#        HHTHSIncrReduct_TP.h
#        HHTGeneralizedExplicit_TP.h
#        NewmarkHSIncrLimit.h
#        NewmarkHSIncrReduct.h
#        NewmarkHSFixedNumIter.h
# )
add_subdirectory(Dynamic)
add_subdirectory(Static)
add_subdirectory(Hybrid)
target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})

