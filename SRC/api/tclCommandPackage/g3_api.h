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

Domain* G3_getDomain(Tcl_Interp*);
TclSafeBuilder *G3_getSafeBuilder(Tcl_Interp *);

UniaxialMaterial* G3_getUniaxialMaterialInstance(Tcl_Interp*, int);
int G3_addUniaxialMaterial(Tcl_Interp*, UniaxialMaterial*);


StaticAnalysis* G3_getStaticAnalysis(Tcl_Interp*);
int G3_setStaticAnalysis(Tcl_Interp*, StaticAnalysis*);
int G3_delStaticAnalysis(Tcl_Interp*);
DirectIntegrationAnalysis* G3_getTransientAnalysis(Tcl_Interp*);
int G3_setTransientAnalysis(Tcl_Interp*, DirectIntegrationAnalysis*);


StaticIntegrator* G3_getStaticIntegrator(Tcl_Interp*);
int G3_setStaticIntegrator(Tcl_Interp*, StaticIntegrator*);


int G3_addTimeSeries(Tcl_Interp *, TimeSeries *);
TimeSeries* G3_getTimeSeries(Tcl_Interp *, int);


#endif // G3_API_H_
