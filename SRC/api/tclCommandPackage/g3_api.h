// Written: cmp

#ifndef G3_API_H_
#define G3_API_H_

#include <tcl.h>

typedef int G3_Tag;
typedef Tcl_Interp G3_Runtime;
#define G3_getDouble Tcl_GetDoubleFromObj

class TclSafeBuilder;
class AnalysisModel;
class EquiSolnAlgo;
class ConstraintHandler;
class DOF_Numberer;
class LinearSOE;
class EigenSOE;
class StaticAnalysis;
class DirectIntegrationAnalysis;
class VariableTimeStepDirectIntegrationAnalysis;
class StaticIntegrator;
class TransientIntegrator;
class ConvergenceTest;

class TimeSeries;

class UniaxialMaterial;
class NDMaterial;
class SectionForceDeformation;
class CrdTransf;
class FrictionModel;
class LimitCurve;
class Domain;
class FE_Datastore;

#ifdef __cplusplus
extern "C" {
#endif

Domain *G3_getDomain(G3_Runtime *);
TclSafeBuilder *G3_getSafeBuilder(G3_Runtime *);

UniaxialMaterial *G3_getUniaxialMaterialInstance(G3_Runtime *, G3_Tag);
int G3_addUniaxialMaterial(G3_Runtime *, UniaxialMaterial *);

StaticAnalysis *G3_getStaticAnalysis(G3_Runtime *);
int G3_setStaticAnalysis(G3_Runtime *, StaticAnalysis *);
int G3_delStaticAnalysis(G3_Runtime *);
DirectIntegrationAnalysis *G3_getTransientAnalysis(G3_Runtime *);
int G3_setTransientAnalysis(G3_Runtime *, DirectIntegrationAnalysis *);

StaticIntegrator *G3_getStaticIntegrator(G3_Runtime *);
int G3_setStaticIntegrator(G3_Runtime *, StaticIntegrator *);

int G3_addTimeSeries(G3_Runtime *, TimeSeries *);
TimeSeries *G3_getTimeSeries(G3_Runtime *, G3_Tag);

#ifdef __cplusplus
}
#endif
#endif // G3_API_H_
