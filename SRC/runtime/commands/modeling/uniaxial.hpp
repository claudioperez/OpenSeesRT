/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Description: This file contains the function invoked when the user invokes
// the uniaxialMaterial command in the interpreter.
//
// Written: fmk, MHS, cmp
// Created: 07/99
//
#include <string>
#include <unordered_map>
#include <runtimeAPI.h>

extern OPS_Routine OPS_SPSW02;           // SAJalali
extern OPS_Routine OPS_TDConcreteEXP;    // ntosic
extern OPS_Routine OPS_TDConcrete;       // ntosic
extern OPS_Routine OPS_TDConcreteMC10;   // ntosic
extern OPS_Routine OPS_TDConcreteMC10NL; // ntosic
extern OPS_Routine OPS_ElasticMaterial;
extern OPS_Routine OPS_ElasticPPMaterial;
extern OPS_Routine OPS_EPPGapMaterial;
// extern void *OPS_ParallelMaterial(G3_Runtime *);
extern OPS_Routine OPS_SeriesMaterial;
extern OPS_Routine OPS_HardeningMaterial;
extern OPS_Routine OPS_HystereticMaterial;
extern OPS_Routine OPS_CableMaterial;
extern OPS_Routine OPS_Bilin;
extern OPS_Routine OPS_Bilin02;
extern OPS_Routine OPS_Steel01;
extern OPS_Routine OPS_FRPConfinedConcrete02;
extern OPS_Routine OPS_Steel02;
extern OPS_Routine OPS_SteelFractureDI; // galvisf
extern OPS_Routine OPS_Steel02Fatigue;
extern OPS_Routine OPS_RambergOsgoodSteel;
extern OPS_Routine OPS_ReinforcingSteel;
extern OPS_Routine OPS_Concrete01;
extern OPS_Routine OPS_Concrete02;
extern OPS_Routine OPS_Concrete02IS;
extern OPS_Routine OPS_PinchingLimitStateMaterial;
extern OPS_Routine OPS_SAWSMaterial;
extern OPS_Routine OPS_ConcreteZ01Material;
extern OPS_Routine OPS_ConcreteL01Material;
extern OPS_Routine OPS_SteelZ01Material;
extern OPS_Routine OPS_TendonL01Material;
extern OPS_Routine OPS_ConfinedConcrete01Material;
extern OPS_Routine OPS_ElasticBilin;
extern OPS_Routine OPS_MinMaxMaterial;
extern OPS_Routine OPS_SimpleFractureMaterial;
extern OPS_Routine OPS_HoehlerStanton;
extern OPS_Routine OPS_InitStrainMaterial;
extern OPS_Routine OPS_InitStressMaterial;
extern OPS_Routine OPS_pyUCLA;
extern OPS_Routine OPS_Maxwell;
extern OPS_Routine OPS_ViscousDamper;
extern OPS_Routine OPS_DamperMaterial;
extern OPS_Routine OPS_BilinearOilDamper;
extern OPS_Routine OPS_Cast;
extern OPS_Routine OPS_Dodd_Restrepo;
extern OPS_Routine OPS_DoddRestr;
extern OPS_Routine OPS_ElasticMultiLinear;
extern OPS_Routine OPS_ImpactMaterial;
extern OPS_Routine OPS_SteelBRB;
extern OPS_Routine OPS_MultiLinear;
extern OPS_Routine OPS_HookGap;
extern OPS_Routine OPS_HyperbolicGapMaterial;
extern OPS_Routine OPS_FRPConfinedConcrete;
extern OPS_Routine OPS_FRPConfinedConcrete02;
extern OPS_Routine OPS_UVCuniaxial;
extern OPS_Routine OPS_Steel01Thermal;
extern OPS_Routine OPS_Steel02Thermal;
extern OPS_Routine OPS_Concrete02Thermal;
extern OPS_Routine OPS_StainlessECThermal;
extern OPS_Routine OPS_SteelECThermal;
extern OPS_Routine OPS_ConcreteECThermal;
extern OPS_Routine OPS_ElasticMaterialThermal;
extern OPS_Routine OPS_ASD_SMA_3K;
extern OPS_Routine OPS_BWBN;
extern OPS_Routine OPS_IMKPeakOriented;
extern OPS_Routine OPS_IMKBilin;
extern OPS_Routine OPS_IMKPinching;
extern OPS_Routine OPS_ModIMKPeakOriented;
extern OPS_Routine OPS_ModIMKPeakOriented02;
extern OPS_Routine OPS_ModIMKPinching;
extern void *OPS_ConcretewBeta(void);
extern OPS_Routine OPS_ConcreteD;
extern OPS_Routine OPS_PinchingLimitState;
extern OPS_Routine OPS_OriginCentered;
extern OPS_Routine OPS_Steel2;
extern OPS_Routine OPS_ConcreteSakaiKawashima;
extern OPS_Routine OPS_ResilienceMaterialHR;
extern OPS_Routine OPS_CFSSSWP;
extern OPS_Routine OPS_CFSWSWP;
extern OPS_Routine OPS_ResilienceLow;
extern OPS_Routine OPS_ViscousMaterial;
extern OPS_Routine OPS_SteelMPF;
extern OPS_Routine OPS_ConcreteCM;
extern OPS_Routine OPS_Bond_SP01;
extern OPS_Routine OPS_Steel4;
extern OPS_Routine OPS_PySimple3;
extern OPS_Routine OPS_BoucWenOriginal;
extern OPS_Routine OPS_GNGMaterial;
extern OPS_Routine OPS_OOHystereticMaterial;
extern OPS_Routine OPS_ElasticPowerFunc;
extern OPS_Routine OPS_UVCuniaxial;
extern OPS_Routine OPS_DegradingPinchedBW;
extern OPS_Routine OPS_SLModel;
extern OPS_Routine OPS_SMAMaterial;
extern OPS_Routine OPS_HystereticPoly;
extern OPS_Routine OPS_Masonry;
extern OPS_Routine OPS_Trilinwp;
extern OPS_Routine OPS_Trilinwp2;
extern OPS_Routine OPS_Masonryt;

const std::unordered_map<std::string, OPS_Routine*> uniaxial_rt_table 
{
// Standard
    {"Elastic",                OPS_ElasticMaterial           },

    {"Concrete01",             OPS_Concrete01                },
    {"Concrete02",             OPS_Concrete02                },

// Composites
    {"MinMaxMaterial",         OPS_MinMaxMaterial            },
    {"MinMax",                 OPS_MinMaxMaterial            },

    // {"Parallel",               OPS_ParallelMaterial          },

    {"Series",                 OPS_SeriesMaterial            },

// Steels

    {"Steel01",                OPS_Steel01                   },

    {"Steel02",                OPS_Steel02                   },

    {"SteelBRB",               OPS_SteelBRB                  },

    {"SteelFractureDI",        OPS_SteelFractureDI           },

    {"Steel02Fatigue",         OPS_Steel02Fatigue            },

    {"Steel4",                 OPS_Steel4                    },

    {"Dodd_Restrepo",          OPS_Dodd_Restrepo             },
    {"DoddRestrepo" ,          OPS_Dodd_Restrepo             },
    {"Restrepo",               OPS_Dodd_Restrepo             },

#if !defined(_NO_NEW_RESTREPO)
    {"DoddRestr",              OPS_DoddRestr                 },
#endif


// Piles
    {"PySimple3",              OPS_PySimple3                 },

// Other
    {"ElasticBilin",           OPS_ElasticBilin              },
    {"ElasticBilinear",        OPS_ElasticBilin              },

    {"ImpactMaterial",         OPS_ImpactMaterial            },
    {"Impact",                 OPS_ImpactMaterial            },

    {"UVCuniaxial",            OPS_UVCuniaxial               },
    {"GNG",                    OPS_GNGMaterial               },

    {"SimpleFractureMaterial", OPS_SimpleFractureMaterial    },
    {"SimpleFracture",         OPS_SimpleFractureMaterial    },

    {"Maxwell",                OPS_Maxwell                   },
    {"MaxwellMaterial",        OPS_Maxwell                   },

    {"ViscousDamper",          OPS_ViscousDamper             },

    {"DamperMaterial",         OPS_DamperMaterial            },

// Concretes
    {"Concrete02IS",           OPS_Concrete02IS              },
    {"ConcreteCM",             OPS_ConcreteCM                },
    {"ConfinedConcrete01",     OPS_ConfinedConcrete01Material},
    {"ConfinedConcrete",       OPS_ConfinedConcrete01Material},

    {"BilinearOilDamper",      OPS_BilinearOilDamper         },

    {"Cast",                   OPS_Cast                      },
    {"CastFuse",               OPS_Cast                      },

    {"ElasticMultiLinear",     OPS_ElasticMultiLinear        },
    {"ElasticPowerFunc",       OPS_ElasticPowerFunc          },

/* 
    {"HoehlerStanton",         OPS_HoehlerStanton            },
*/  

    {"SLModel",                OPS_SLModel                   },

    {"RambergOsgood",          OPS_RambergOsgoodSteel        },
    {"RambergOsgoodSteel",     OPS_RambergOsgoodSteel        },

//  {"ReinforcingSteel",       OPS_ReinforcingSteel          },

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

    {"Hardening",              OPS_HardeningMaterial         },
    {"Hardening2",             OPS_HardeningMaterial         },

    {"BilinMaterial",          OPS_Bilin                     },
    {"Bilin",                  OPS_Bilin                     },
    
    {"Hysteretic",             OPS_HystereticMaterial        },

    {"ElasticPPGap",           OPS_EPPGapMaterial            },


    {"OOHysteretic",           OPS_OOHystereticMaterial      },

    {"Viscous",                OPS_ViscousMaterial           },

    {"SAWSMaterial",           OPS_SAWSMaterial              },
    {"SAWS",                   OPS_SAWSMaterial              },

    {"ConcreteZ01Material",    OPS_ConcreteZ01Material       },
    {"ConcreteZ01",            OPS_ConcreteZ01Material       },

    {"ConcreteL01Material",    OPS_ConcreteL01Material       },
    {"ConcreteL01",            OPS_ConcreteL01Material       },

    {"SteelZ01Material",       OPS_SteelZ01Material          },
    {"SteelZ01",               OPS_SteelZ01Material          },

    {"TendonL01Material",      OPS_TendonL01Material         },
    {"TendonL01",              OPS_TendonL01Material         },

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

typedef UniaxialMaterial* (G3_TclUniaxialCommand)(G3_Runtime*, int, TCL_Char ** const);
G3_TclUniaxialCommand TclCommand_KikuchiAikenHDR;
G3_TclUniaxialCommand TclCommand_KikuchiAikenLRB;
G3_TclUniaxialCommand G3Parse_newFedeasUniaxialDamage;
G3_TclUniaxialCommand G3Parse_newUniaxialConcrete04;
G3_TclUniaxialCommand G3Parse_newUniaxialConcrete06;
G3_TclUniaxialCommand G3Parse_newUniaxialConcrete07;
G3_TclUniaxialCommand TclCommand_ReinforcingSteel;
G3_TclUniaxialCommand G3Parse_newParallelMaterial;
G3_TclUniaxialCommand G3Parse_newUniaxialBoucWen;
// G3_TclUniaxialCommand TclCommand_AxialSp;
// G3_TclUniaxialCommand TclCommand_AxialSpHD;

template <void*(*fn)(G3_Runtime*)> static void*
G3_(G3_Runtime* rt, int argc, G3_Char** const)
{
  return fn(rt);
}


std::unordered_map<std::string, G3_TclUniaxialCommand *> uniaxial_tcl_table = {

    {"FedeasUniaxialDamage", G3Parse_newFedeasUniaxialDamage  },
    {"KikuchiAikenHDR",      TclCommand_KikuchiAikenHDR       },
    {"KikuchiAikenLRB",      TclCommand_KikuchiAikenLRB       },
    /*
    {"AxialSp",             TclCommand_AxialSp               },
    {"AxialSpHD",           TclCommand_AxialSpHD             },
    */
    {"Concrete04",          G3Parse_newUniaxialConcrete04 },
    {"Concrete06",          G3Parse_newUniaxialConcrete06 },
    {"Concrete07",          G3Parse_newUniaxialConcrete07 },
    {"ReinforcingSteel",    TclCommand_ReinforcingSteel      }, 
    {"Parallel",            G3Parse_newParallelMaterial      },
    {"BoucWen",             G3Parse_newUniaxialBoucWen       },
/*
    {"Elastic",             G3_<OPS_ElasticMaterial>         },

    {"Steel01",             G3_<OPS_Steel01>                 },

    {"Steel02",             G3_<OPS_Steel02>                 }
*/
};


typedef UniaxialMaterial*(G3_TclUniaxialPackage)(ClientData, Tcl_Interp *, int, TCL_Char ** const);
G3_TclUniaxialPackage TclBasicBuilder_addFedeasMaterial;
G3_TclUniaxialPackage TclBasicBuilder_addSnapMaterial;
G3_TclUniaxialPackage TclBasicBuilder_addDrainMaterial;
std::unordered_map<std::string, G3_TclUniaxialPackage *> tcl_uniaxial_package_table {

  {"DRAIN",              TclBasicBuilder_addDrainMaterial },

  {"SNAP",               TclBasicBuilder_addSnapMaterial  },
  {"snap",               TclBasicBuilder_addSnapMaterial  },

// #if defined(_STEEL2) || defined(OPSDEF_UNIAXIAL_FEDEAS)
//{"FEDEAS",             TclBasicBuilder_addFedeasMaterial},
// #endif
};

