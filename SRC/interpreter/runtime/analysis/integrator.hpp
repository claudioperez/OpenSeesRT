#include <string>
#include <unordered_map>
#include <g3_api.h>
#include <InputAPI.h>
#include <runtimeAPI.h>

#include <Domain.h>
#include <Node.h>

// integrators
#include <LoadControl.h>
#include <StagedLoadControl.h>
#include <ArcLength1.h>
#include <DisplacementControl.h>

#include <Newmark.h>
#include <StagedNewmark.h>
#include <TRBDF2.h>
#include <TRBDF3.h>
#include <Houbolt.h>
#include <ParkLMS3.h>
#include <BackwardEuler.h>

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char **argv, Domain *domain);

StaticIntegrator*
G3Parse_newHSIntegrator(G3_Builder *, int, const char **);
StaticIntegrator*
G3Parse_newLoadControl(G3_Builder *, int argc, const char *argv[]);
StaticIntegrator*
G3Parse_newEQPathIntegrator(G3_Builder *, int argc, const char *argv[]);
StaticIntegrator*
G3Parse_newArcLengthIntegrator(G3_Builder *, int argc, const char *argv[]);
StaticIntegrator*
G3Parse_newStagedLoadControlIntegrator(G3_Builder*, int, TCL_Char **);
StaticIntegrator*
G3Parse_newMinUnbalDispNormIntegrator(G3_Runtime*, int, G3_Char **);
StaticIntegrator*
G3Parse_newDisplacementControlIntegrator(G3_Builder *, int, G3_Char**);
StaticIntegrator*
G3Parse_newStaticIntegrator(G3_Builder *, int, TCL_Char **);

TransientIntegrator*
G3Parse_newNewmark1Integrator(G3_Builder *, int, TCL_Char **);
TransientIntegrator*
G3Parse_newNewmarkIntegrator(G3_Runtime*, int, G3_Char**);
TransientIntegrator*
G3Parse_newTransientIntegrator(G3_Builder *, int, TCL_Char **);

std::unordered_map<std::string, StaticIntegrator* (*)(G3_Runtime*, int, G3_Char**)> 
StaticIntegratorLibrary = {
  {"LoadControl",                  G3Parse_newLoadControl},
  {"StagedLoadControl",            G3Parse_newStagedLoadControlIntegrator},
  {"EQPath",                       G3Parse_newEQPathIntegrator},
  {"ArcLength",                    G3Parse_newArcLengthIntegrator},
  {"MinUnbalDispNorm",             G3Parse_newMinUnbalDispNormIntegrator},
  {"DisplacementControl",          G3Parse_newDisplacementControlIntegrator},
};

typedef TransientIntegrator* (*TransientConstructor)(G3_Runtime*, int, G3_Char**);
std::unordered_map<std::string, TransientIntegrator* (*)(G3_Runtime*, int, G3_Char**)> 
TransientIntegratorLibrary = {
  {"TRBDF2",         (TransientConstructor)[](G3_Runtime*, int, G3_Char**)->TransientIntegrator*{return new TRBDF2();}},
  {"Bathe",          (TransientConstructor)[](G3_Runtime*, int, G3_Char**)->TransientIntegrator*{return new TRBDF2();}},
  {"TRBDF3",         (TransientConstructor)[](G3_Runtime*, int, G3_Char**)->TransientIntegrator*{return new TRBDF3();}},
  {"Bathe3",         (TransientConstructor)[](G3_Runtime*, int, G3_Char**)->TransientIntegrator*{return new TRBDF3();}},
  {"Houbolt",        (TransientConstructor)[](G3_Runtime*, int, G3_Char**)->TransientIntegrator*{return new Houbolt();}},
  {"ParkLMS3",       (TransientConstructor)[](G3_Runtime*, int, G3_Char**)->TransientIntegrator*{return new ParkLMS3();}},

  {"Newmark",        G3Parse_newNewmarkIntegrator},
};

