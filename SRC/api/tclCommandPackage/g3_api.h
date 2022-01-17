// Written: cmp

#ifndef G3_API_H_
#define G3_API_H_

#define G3_MAX_NUM_DOFS 1000000000000
#define G3_NUM_DOF_BUFFER 20
#include <tcl.h>
// #include <g3_io.h>

#ifndef OPS_Export
# define OPS_Export
#endif

#define OPS_GetDomain() G3_getDomain(rt)
#define OPS_GetUniaxialMaterial(tag) G3_getUniaxialMaterialInstance(rt, (tag))

typedef int G3_Tag;
// typedef Tcl_Interp G3_Runtime;
// #define G3_Runtime Tcl_Interp

#define G3_getDouble Tcl_GetDoubleFromObj
class G3_Runtime;
class ModelBuilder;
class TclBuilder;
class TclSafeBuilder;
class TclBasicBuilder;
class AnalysisModel;
class EquiSolnAlgo;
class ConstraintHandler;
class DOF_Numberer;
class LinearSOE;
class EigenSOE;
class AnalysisModel;
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

typedef int (G3_TclElementCommand)(ClientData, Tcl_Interp*, int, const char**, Domain*, TclBasicBuilder*);

#ifdef __cplusplus
extern "C" {
#endif

// Runtime
G3_Runtime *G3_getRuntime(Tcl_Interp *);
Tcl_Interp *G3_getInterpreter(G3_Runtime*);

// Domain
int G3_setDomain(G3_Runtime*, Domain*);
Domain *G3_getDomain(G3_Runtime *);

TclSafeBuilder *G3_getSafeBuilder(G3_Runtime *);
TclBuilder     *G3_getModelBuilder(G3_Runtime *);
int             G3_setModelBuilder(G3_Runtime *, TclBuilder*);
bool G3_modelIsBuilt(G3_Runtime *);
int G3_getNDM(G3_Runtime *);
int G3_getNDF(G3_Runtime *);

// Materials
UniaxialMaterial *G3_getUniaxialMaterialInstance(G3_Runtime *, G3_Tag);
int G3_addUniaxialMaterial(G3_Runtime *, UniaxialMaterial *);


// Analysis
AnalysisModel *G3_getAnalysisModel(G3_Runtime *);
int G3_setAnalysisModel(G3_Runtime *, AnalysisModel *);

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
#include <elementAPI.h>
#endif // G3_API_H_
