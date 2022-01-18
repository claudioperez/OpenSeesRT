add_library(g3_uniaxial_io OBJECT)
target_include_directories(g3_uniaxial_io PRIVATE ${PROJECT_SOURCE_DIR})

target_sources(g3_uniaxial_io
    PRIVATE
        IO_ASD_SMA_3K.cpp
        IO_AxialSp.cpp
        IO_AxialSpHD.cpp
        IO_BackboneMaterial.cpp
        IO_BarSlipMaterial.cpp
        # IO_Bilin02.cpp
        # IO_Bilin.cpp
        # IO_BilinearOilDamper.cpp
        IO_Bond_SP01.cpp
        IO_BoucWenMaterial.cpp
        IO_BoucWenOriginal.cpp
        IO_BraceMaterial.cpp
        IO_BWBN.cpp
        IO_CableMaterial.cpp
        # IO_Cast.cpp
        # IO_CFSSSWP.cpp
        # IO_CFSWSWP.cpp
        IO_Concrete01A.cpp
        IO_Concrete01.cpp
        IO_Concrete01WithSITC.cpp
        IO_Concrete02.cpp
        IO_Concrete02IS.cpp
        IO_Concrete02Thermal.cpp
        IO_Concrete04.cpp
        IO_Concrete05.cpp
        IO_Concrete06.cpp
        IO_Concrete07.cpp
        IO_ConcreteCM.cpp
        # IO_ConcreteD.cpp
        IO_ConcreteECThermal.cpp
        IO_ConcreteSakaiKawashima.cpp
        IO_ConcretewBeta.cpp
        # IO_ConfinedConcrete01.cpp
        IO_ContinuumUniaxial.cpp
        IO_CubicSpline.cpp
        IO_DamperMaterial.cpp
        IO_DegradingPinchedBW.cpp
        # IO_DoddRestr.cpp
        # IO_Dodd_Restrepo.cpp
        IO_DrainMaterial.cpp
        IO_ECC01.cpp
        IO_Elastic2Material.cpp
        IO_ElasticBDMaterial.cpp
        IO_ElasticBilin.cpp
        IO_ElasticMaterial.cpp
        IO_ElasticMaterialThermal.cpp
        IO_ElasticMultiLinear.cpp
        IO_ElasticPowerFunc.cpp
        IO_ElasticPPMaterial.cpp
        IO_ENTMaterial.cpp
        IO_EPPGapMaterial.cpp
        # IO_ExternalUniaxialMaterial.cpp
        IO_FatigueMaterial.cpp
        IO_FedeasMaterial.cpp
        # IO_FRPConfinedConcrete02.cpp
        # IO_FRPConfinedConcrete.cpp
        # IO_GNGMaterial.cpp
        IO_HardeningMaterial2.cpp
        IO_HardeningMaterial.cpp
        IO_HookGap.cpp
        IO_HyperbolicGapMaterial.cpp
        IO_HystereticMaterial.cpp
        IO_HystereticPoly.cpp
        # IO_IMKBilin.cpp
        # IO_IMKPeakOriented.cpp
        # IO_IMKPinching.cpp
        IO_ImpactMaterial.cpp
        IO_InitStrainMaterial.cpp
        IO_InitStressMaterial.cpp
        IO_KikuchiAikenHDR.cpp
        IO_KikuchiAikenLRB.cpp
        # IO_Masonry2.cpp
        # IO_Masonry.cpp
        IO_Masonryt.cpp
        IO_MaterialState.cpp
        # IO_Maxwell.cpp
        IO_MinMaxMaterial.cpp
        # IO_ModIMKPeakOriented02.cpp
        # IO_ModIMKPeakOriented.cpp
        # IO_ModIMKPinching02.cpp
        # IO_ModIMKPinching.cpp
        IO_MultiLinear.cpp
        IO_MultiplierMaterial.cpp
        IO_Neoprene.cpp
        IO_NewUniaxialMaterial.cpp
        # IO_OOHystereticMaterial.cpp
        IO_OriginCentered.cpp
        IO_ParallelMaterial.cpp
        IO_PathIndependentMaterial.cpp
        IO_PenaltyMaterial.cpp
        IO_Pinching4Material.cpp
        IO_PrestressedSteelMaterial.cpp
        # IO_pyUCLA.cpp
        # IO_RambergOsgoodSteel.cpp
        IO_ReinforcingSteel.cpp
        IO_ResilienceLow.cpp
        IO_ResilienceMaterialHR.cpp
        # IO_SAWSMaterial.cpp
        IO_SecantConcrete.cpp
        IO_SelfCenteringMaterial.cpp
        IO_SeriesMaterial.cpp
        IO_ShearPanelMaterial.cpp
        IO_SimpleFractureMaterial.cpp
        # IO_SLModel.cpp
        IO_SMAMaterial.cpp
        IO_SmoothPSConcrete.cpp
        # IO_SPSW02.cpp
        IO_StainlessECThermal.cpp
        IO_Steel01.cpp
        # IO_Steel01Thermal.cpp
        IO_Steel02.cpp
        IO_Steel02Fatigue.cpp
        IO_Steel02Thermal.cpp
        IO_Steel03.cpp
        IO_Steel2.cpp
        IO_Steel4.cpp
        IO_SteelBRB.cpp
        IO_SteelECThermal.cpp
        IO_SteelFractureDI.cpp
        IO_SteelMP.cpp
        IO_SteelMPF.cpp
        # IO_TDConcrete.cpp
        # IO_TDConcreteEXP.cpp
        IO_TDConcreteMC10.cpp
        IO_TDConcreteMC10NL.cpp
        IO_TensionOnlyMaterial.cpp
        IO_Trilinwp2.cpp
        IO_Trilinwp.cpp
        # IO_Trilinwpd.cpp
        IO_TriMatrix.cpp
        IO_UniaxialJ2Plasticity.cpp
        IO_UniaxialMaterial.cpp
        # IO_UVCuniaxial.cpp
        # IO_ViscousDamper.cpp
        # IO_ViscousMaterial.cpp
        IO_WrapperUniaxialMaterial.cpp
)

#add_subdirectory(./backbone/)
#add_subdirectory(./drain/)
#add_subdirectory(./fedeas/)
#add_subdirectory(./limitState/)
#add_subdirectory(./PY/)
#add_subdirectory(./snap/)
#add_subdirectory(./stiffness/)
#add_subdirectory(./strength/)
#add_subdirectory(./unloading/)

