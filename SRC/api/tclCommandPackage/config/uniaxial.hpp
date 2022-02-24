// Written: fmk, MHS, cmp
//
// Description: This file contains the function invoked when the user invokes
// the uniaxialMaterial command in the interpreter.
//
#include <string>
#include <unordered_map>
extern void *OPS_SPSW02(G3_Runtime *);           // SAJalali
extern void *OPS_TDConcreteEXP(G3_Runtime *);    // ntosic
extern void *OPS_TDConcrete(G3_Runtime *);       // ntosic
extern void *OPS_TDConcreteMC10(G3_Runtime *);   // ntosic
extern void *OPS_TDConcreteMC10NL(G3_Runtime *); // ntosic
extern void *OPS_ElasticMaterial(G3_Runtime *);
extern void *OPS_ElasticPPMaterial(G3_Runtime *);
extern void *OPS_EPPGapMaterial(G3_Runtime *);
extern void *OPS_ParallelMaterial(G3_Runtime *);
extern void *OPS_SeriesMaterial(G3_Runtime *);
extern void *OPS_HardeningMaterial(G3_Runtime *);
extern void *OPS_HystereticMaterial(G3_Runtime *);
extern void *OPS_CableMaterial(G3_Runtime *);
extern void *OPS_Bilin(G3_Runtime *);
extern void *OPS_Bilin02(G3_Runtime *);
extern void *OPS_Steel01(G3_Runtime *);
extern void *OPS_FRPConfinedConcrete02(G3_Runtime *);
extern void *OPS_Steel02(G3_Runtime *);
extern void *OPS_SteelFractureDI(G3_Runtime *); // galvisf
extern void *OPS_Steel02Fatigue(G3_Runtime *);
extern void *OPS_RambergOsgoodSteel(G3_Runtime *);
extern void *OPS_ReinforcingSteel(G3_Runtime *);
extern void *OPS_Concrete01(G3_Runtime *);
extern void *OPS_Concrete02(G3_Runtime *);
extern void *OPS_Concrete02IS(G3_Runtime *);
extern void *OPS_PinchingLimitStateMaterial(G3_Runtime *);
extern void *OPS_SAWSMaterial(G3_Runtime *);
extern void *OPS_ConcreteZ01Material(G3_Runtime *);
extern void *OPS_ConcreteL01Material(G3_Runtime *);
extern void *OPS_SteelZ01Material(G3_Runtime *);
extern void *OPS_TendonL01Material(G3_Runtime *);
extern void *OPS_ConfinedConcrete01Material(G3_Runtime *);
extern void *OPS_ElasticBilin(G3_Runtime *);
extern void *OPS_MinMaxMaterial(G3_Runtime *);
extern void *OPS_SimpleFractureMaterial(G3_Runtime *);
extern void *OPS_HoehlerStanton(G3_Runtime *);
extern void *OPS_InitStrainMaterial(G3_Runtime *);
extern void *OPS_InitStressMaterial(G3_Runtime *);
extern void *OPS_pyUCLA(G3_Runtime *);
extern void *OPS_Maxwell(G3_Runtime *);
extern void *OPS_ViscousDamper(G3_Runtime *);
extern void *OPS_DamperMaterial(G3_Runtime *);
extern void *OPS_BilinearOilDamper(G3_Runtime *);
extern void *OPS_Cast(G3_Runtime *);
extern void *OPS_Dodd_Restrepo(G3_Runtime *);
extern void *OPS_DoddRestr(G3_Runtime *);
extern void *OPS_ElasticMultiLinear(G3_Runtime *);
extern void *OPS_ImpactMaterial(G3_Runtime *);
extern void *OPS_SteelBRB(G3_Runtime *);
extern void *OPS_MultiLinear(G3_Runtime *);
extern void *OPS_HookGap(G3_Runtime *);
extern void *OPS_HyperbolicGapMaterial(G3_Runtime *);
extern void *OPS_FRPConfinedConcrete(G3_Runtime *);
extern void *OPS_FRPConfinedConcrete02(G3_Runtime *);
extern void *OPS_UVCuniaxial(G3_Runtime *);
extern void *OPS_Steel01Thermal(G3_Runtime *);
extern void *OPS_Steel02Thermal(G3_Runtime *);
extern void *OPS_Concrete02Thermal(G3_Runtime *);
extern void *OPS_StainlessECThermal(G3_Runtime *);
extern void *OPS_SteelECThermal(G3_Runtime *);
extern void *OPS_ConcreteECThermal(G3_Runtime *);
extern void *OPS_ElasticMaterialThermal(G3_Runtime *);
extern void *OPS_ASD_SMA_3K(G3_Runtime *);
extern void *OPS_BWBN(G3_Runtime *);
extern void *OPS_IMKPeakOriented(G3_Runtime *);
extern void *OPS_IMKBilin(G3_Runtime *);
extern void *OPS_IMKPinching(G3_Runtime *);
extern void *OPS_ModIMKPeakOriented(G3_Runtime *);
extern void *OPS_ModIMKPeakOriented02(G3_Runtime *);
extern void *OPS_ModIMKPinching(G3_Runtime *);
extern void *OPS_ConcretewBeta(void);
extern void *OPS_ConcreteD(G3_Runtime *);
extern void *OPS_PinchingLimitState(G3_Runtime *);
extern void *OPS_OriginCentered(G3_Runtime *);
extern void *OPS_Steel2(G3_Runtime *);
extern void *OPS_ConcreteSakaiKawashima(G3_Runtime *);
extern void *OPS_ResilienceMaterialHR(G3_Runtime *);
extern void *OPS_CFSSSWP(G3_Runtime *);
extern void *OPS_CFSWSWP(G3_Runtime *);
extern void *OPS_ResilienceLow(G3_Runtime *);
extern void *OPS_ViscousMaterial(G3_Runtime *);
extern void *OPS_SteelMPF(G3_Runtime *);
extern void *OPS_ConcreteCM(G3_Runtime *);
extern void *OPS_Bond_SP01(G3_Runtime *);
extern void *OPS_Steel4(G3_Runtime *);
extern void *OPS_PySimple3(G3_Runtime *);
extern void *OPS_BoucWenOriginal(G3_Runtime *);
extern void *OPS_GNGMaterial(G3_Runtime *);
extern void *OPS_OOHystereticMaterial(G3_Runtime *);
extern void *OPS_ElasticPowerFunc(G3_Runtime *);
extern void *OPS_UVCuniaxial(G3_Runtime *);
extern void *OPS_DegradingPinchedBW(G3_Runtime *);
extern void *OPS_SLModel(G3_Runtime *);
extern void *OPS_SMAMaterial(G3_Runtime *);
extern void *OPS_HystereticPoly(G3_Runtime *);
extern void *OPS_Masonry(G3_Runtime *);
extern void *OPS_Trilinwp(G3_Runtime *);
extern void *OPS_Trilinwp2(G3_Runtime *);
extern void *OPS_Masonryt(G3_Runtime *);

typedef void* (G3_RuntimeUniaxialCommand)(G3_Runtime*);

const std::unordered_map<std::string, G3_RuntimeUniaxialCommand*> uniaxial_rt_table 
{

    {"Elastic",                OPS_ElasticMaterial           },

    {"Steel01",                OPS_Steel01                   },

    {"Steel02",                OPS_Steel02                   },

    {"Parallel",               OPS_ParallelMaterial          },

    {"Series",                 OPS_SeriesMaterial            },

    {"SteelFractureDI",        OPS_SteelFractureDI           },

    {"Steel02Fatigue",         OPS_Steel02Fatigue            },
    {"Steel4",                 OPS_Steel4                    },
    {"UVCuniaxial",            OPS_UVCuniaxial               },
    {"GNG",                    OPS_GNGMaterial               },
    {"PySimple3",              OPS_PySimple3                 },
    {"Concrete01",             OPS_Concrete01                },
    {"Concrete02",             OPS_Concrete02                },
    {"Concrete02IS",           OPS_Concrete02IS              },

    {"ElasticBilin",           OPS_ElasticBilin              },
    {"ElasticBilinear",        OPS_ElasticBilin              },

    {"ImpactMaterial",         OPS_ImpactMaterial            },
    {"Impact",                 OPS_ImpactMaterial            },

    {"SteelBRB",               OPS_SteelBRB                  },

    {"MinMaxMaterial",         OPS_MinMaxMaterial            },
    {"MinMax",                 OPS_MinMaxMaterial            },

    {"SimpleFractureMaterial", OPS_SimpleFractureMaterial    },
    {"SimpleFracture",         OPS_SimpleFractureMaterial    },

    {"Maxwell",                OPS_Maxwell                   },
    {"MaxwellMaterial",        OPS_Maxwell                   },

    {"ViscousDamper",          OPS_ViscousDamper             },

    {"DamperMaterial",         OPS_DamperMaterial            },

    {"BilinearOilDamper",      OPS_BilinearOilDamper         },

    {"Cast",                   OPS_Cast                      },
    {"CastFuse",               OPS_Cast                      },

    {"Dodd_Restrepo",          OPS_Dodd_Restrepo             },
    {"DoddRestrepo" ,          OPS_Dodd_Restrepo             },
    {"Restrepo",               OPS_Dodd_Restrepo             },

#if !defined(_NO_NEW_RESTREPO)
    {"DoddRestr",              OPS_DoddRestr                 },
#endif

    {"ElasticMultiLinear",     OPS_ElasticMultiLinear        },
    {"ElasticPowerFunc",       OPS_ElasticPowerFunc          },

/* 
    {"HoehlerStanton",         OPS_HoehlerStanton            },
*/  

    {"SLModel",                OPS_SLModel                   },

    {"RambergOsgood",          OPS_RambergOsgoodSteel        },
    {"RambergOsgoodSteel",     OPS_RambergOsgoodSteel        },

    {"ReinforcingSteel",       OPS_ReinforcingSteel          },

    {"Steel2",                 OPS_Steel2                    },

    {"OriginCentered",         OPS_OriginCentered            },

    {"HookGap",                OPS_HookGap                   },

    {"HyperbolicGapMaterial",  OPS_HyperbolicGapMaterial     },

    {"FRPConfinedConcrete02",  OPS_FRPConfinedConcrete02     },

    {"PinchingLimitState",     OPS_PinchingLimitState        },

    {"InitStrainMaterial",     OPS_InitStrainMaterial        },
    {"InitStrain",             OPS_InitStrainMaterial        },

    {"InitStressMaterial",     OPS_InitStressMaterial        },
    {"InitStress",             OPS_InitStressMaterial        },

    {"pyUCLA",                 OPS_pyUCLA                    },
    {"PYUCLA",                 OPS_pyUCLA                    },

    {"MultiLinear",            OPS_MultiLinear               },

    {"BWBN",                   OPS_BWBN                      },

    {"DegradingPinchedBW",     OPS_DegradingPinchedBW        },

    {"IMKBilin",               OPS_IMKBilin                  },

    {"IMKPeakOriented",        OPS_IMKPeakOriented           },

    {"IMKPinching",            OPS_IMKPinching               },

    {"ModIMKPeakOriented",     OPS_ModIMKPeakOriented        },

    {"ModIMKPeakOriented02",   OPS_ModIMKPeakOriented02      },

    {"Bilin02",                OPS_Bilin02                   },

    {"BoucWenOriginal",        OPS_BoucWenOriginal           },

// Thermal
    {"Steel01Thermal",         OPS_Steel01Thermal            },

    {"Steel02Thermal",         OPS_Steel02Thermal            },

    {"SteelECThermal",         OPS_SteelECThermal            },

    {"StainlessECThermal",     OPS_StainlessECThermal        },

    {"ElasticThermal",         OPS_ElasticMaterialThermal    },

    {"ConcreteECThermal",      OPS_ConcreteECThermal         },

    {"Concrete02Thermal",      OPS_Concrete02Thermal         },



#if 0
    {"ConcretewBeta",          OPS_ConcretewBeta             },
#endif
    {"ConcreteD",              OPS_ConcreteD                 },

    {"ConcreteSakaiKawashima", OPS_ConcreteSakaiKawashima    },


    {"SteelMPF",               OPS_SteelMPF                  },

    {"ConcreteCM",             OPS_ConcreteCM                },

    {"ResilienceLow",          OPS_ResilienceLow             },

    {"ResilienceMaterialHR",   OPS_ResilienceMaterialHR      },

    {"CFSWSWP",                OPS_CFSWSWP                   },

    {"CFSSSWP",                OPS_CFSSSWP                   },

    {"FRPConfinedConcrete",    OPS_FRPConfinedConcrete       },

    {"Masonry",                OPS_Masonry                   },

    {"Trilinwp",               OPS_Trilinwp                  },

    {"Trilinwp2",              OPS_Trilinwp2                 },

    {"Masonryt",               OPS_Masonryt                  },

    {"ElasticPP",              OPS_ElasticPPMaterial         },

    {"ElasticPPGap",           OPS_EPPGapMaterial            },

    {"Hardening",              OPS_HardeningMaterial         },
    {"Hardening2",             OPS_HardeningMaterial         },

    {"Hysteretic",             OPS_HystereticMaterial        },

    {"OOHysteretic",           OPS_OOHystereticMaterial      },

    {"Viscous",                OPS_ViscousMaterial           },

    {"SAWSMaterial",           OPS_SAWSMaterial              },
    {"SAWS",                   OPS_SAWSMaterial              },

    {"BilinMaterial",          OPS_Bilin                     },
    {"Bilin",                  OPS_Bilin                     },

    {"ConcreteZ01Material",    OPS_ConcreteZ01Material       },
    {"ConcreteZ01",            OPS_ConcreteZ01Material       },

    {"ConcreteL01Material",    OPS_ConcreteL01Material       },
    {"ConcreteL01",            OPS_ConcreteL01Material       },

    {"SteelZ01Material",       OPS_SteelZ01Material          },
    {"SteelZ01",               OPS_SteelZ01Material          },

    {"TendonL01Material",      OPS_TendonL01Material         },
    {"TendonL01",              OPS_TendonL01Material         },

    {"ConfinedConcrete01",     OPS_ConfinedConcrete01Material},
    {"ConfinedConcrete",       OPS_ConfinedConcrete01Material},

    {"Cable",                  OPS_CableMaterial             },

    {"SMA",                    OPS_SMAMaterial               },

    {"ASD_SMA_3K",             OPS_ASD_SMA_3K                },

    {"HystereticPoly",         OPS_HystereticPoly            },

    {"SPSW02",                 OPS_SPSW02                    },

    {"TDConcreteEXP",          OPS_TDConcreteEXP             },

    {"TDConcrete",             OPS_TDConcrete                },

    {"TDConcreteMC10",         OPS_TDConcreteMC10            },

    {"TDConcreteMC10NL",       OPS_TDConcreteMC10NL          },
};

/*
  {"PlateBearingConnectionThermal",  OPS_PlateBearingConnectionThermal},
  {"PinchingLimitStateMaterial",     OPS_PinchingLimitState           },

*/

typedef int(G3_TclUniaxialCommand)(ClientData, Tcl_Interp *, int, TCL_Char **);
G3_TclUniaxialCommand TclSafeBuilder_addFedeasWrapper;
G3_TclUniaxialCommand TclCommand_KikuchiAikenHDR;
G3_TclUniaxialCommand TclCommand_KikuchiAikenLRB;
G3_TclUniaxialCommand TclCommand_AxialSp;
G3_TclUniaxialCommand TclCommand_AxialSpHD;

std::unordered_map<std::string, G3_TclUniaxialCommand *> uniaxial_tcl_table = {
    {"FedeasDamageWrapper", TclSafeBuilder_addFedeasWrapper  },
    {"KikuchiAikenHDR",     TclCommand_KikuchiAikenHDR       },
    {"KikuchiAikenLRB",     TclCommand_KikuchiAikenLRB       },
    {"AxialSp",             TclCommand_AxialSp               },
    {"AxialSpHD",           TclCommand_AxialSpHD             }
};


typedef UniaxialMaterial*(G3_TclUniaxialPackage)(ClientData, Tcl_Interp *, int, TCL_Char **);
G3_TclUniaxialPackage TclBasicBuilder_addFedeasMaterial;
G3_TclUniaxialPackage TclBasicBuilder_addSnapMaterial;
G3_TclUniaxialPackage TclBasicBuilder_addDrainMaterial;

std::unordered_map<std::string, G3_TclUniaxialPackage *> tcl_uniaxial_package_table {
  {"DRAIN",              TclBasicBuilder_addDrainMaterial },

  {"SNAP",               TclBasicBuilder_addSnapMaterial  },
  {"snap",               TclBasicBuilder_addSnapMaterial  },

// #if defined(_STEEL2) || defined(OPSDEF_UNIAXIAL_FEDEAS)
  {"FEDEAS",             TclBasicBuilder_addFedeasMaterial},
// #endif
};

