/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Written: cmp
//
#ifndef G3_API_H_
#define G3_API_H_
   struct Tcl_Interp;
   typedef void* ClientData;
#  include <api/elementAPI.h>
//
#  ifndef OPS_Export
#    define OPS_Export
#  endif
//
// ANALYSIS RUNTIME
//
#  undef  OPS_GetDomain
// #  define ops_getdomain_
#  define OPS_GetDomain()        G3_getDomain(rt)
#  undef  OPS_GetAnalysisModel
#  define OPS_GetAnalysisModel() G3_getAnalysisModelPtr(rt)
//
// MODELBUILDER
//
// Materials
#  undef  OPS_GetUniaxialMaterial
#  define OPS_GetUniaxialMaterial(tag) G3_getUniaxialMaterialInstance(rt, (tag))
#  define OPS_getUniaxialMaterial(tag) G3_getUniaxialMaterialInstance(rt, (tag))
// Time series
#  define OPS_addTimeSeries(series) G3_addTimeSeries(rt, (series))
#  define OPS_getTimeSeries(tag) G3_getTimeSeries(rt, (tag))
// #  define OPS_removeTimeSeries(tag) G3_removeTimeSeries(rt, (tag))

// Coordinate Transforms
#  undef  OPS_getCrdTransf
#  define OPS_getCrdTransf(tag) G3_getCrdTransf(rt, (tag))
#  undef  OPS_GetCrdTransf
#  define OPS_GetCrdTransf(tag) G3_getCrdTransf(rt, (tag))

// Section
SectionForceDeformation *G3_getSectionForceDeformation(G3_Runtime *, int);
#  undef  OPS_getSectionForceDeformation
#  define OPS_getSectionForceDeformation(tag) G3_getSectionForceDeformation(rt, (tag))
#  undef  OPS_GetSectionForceDeformation
#  define OPS_GetSectionForceDeformation(tag) G3_getSectionForceDeformation(rt, (tag))


typedef int G3_Tag;
#  define G3_Char TCL_Char
#  define G3_getDouble Tcl_GetDoubleFromObj
class G3_Runtime;
class ModelBuilder;
class TclBuilder;
#include <runtime/BasicModelBuilder.h>
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


#ifdef __cplusplus
extern "C" {
#endif
//
// Runtime
//
G3_Runtime *G3_getRuntime(Tcl_Interp *);
Tcl_Interp *G3_getInterpreter(G3_Runtime*);

// Domain
Domain *G3_getDomain(G3_Runtime *);

BasicModelBuilder *G3_getSafeBuilder(G3_Runtime *);
TclBuilder     *G3_getModelBuilder(G3_Runtime *);
int             G3_setModelBuilder(G3_Runtime *, BasicModelBuilder*);
int     G3_getNDM(G3_Runtime *);

// Materials
UniaxialMaterial *G3_getUniaxialMaterialInstance(G3_Runtime *, G3_Tag);
int G3_addUniaxialMaterial(G3_Runtime *, UniaxialMaterial *);
// Coordinate Transforms
CrdTransf *G3_getCrdTransf(G3_Runtime *, G3_Tag);


// Systems and Solvers
LinearSOE **G3_getLinearSoePtr(G3_Runtime* );

//
// RUNTIME
//
// Time Series
int G3_addTimeSeries(G3_Runtime *, TimeSeries *);
TimeSeries *G3_getTimeSeries(G3_Runtime *, G3_Tag);
int G3_removeTimeSeries(G3_Runtime *, G3_Tag);
// Analysis
AnalysisModel **G3_getAnalysisModelPtr(G3_Runtime *);
// StaticAnalysis *G3_getStaticAnalysis(G3_Runtime *);
// int G3_setStaticAnalysis(G3_Runtime *, StaticAnalysis *);
// int G3_delStaticAnalysis(G3_Runtime *);
StaticIntegrator *G3_getStaticIntegrator(G3_Runtime *);
int G3_setStaticIntegrator(G3_Runtime *, StaticIntegrator *);

#ifdef __cplusplus
}
#endif
#endif // G3_API_H_
