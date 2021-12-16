/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

/*
** $Revision: 1.9 $
** $Date: 2010-03-05 22:32:36 $
** $Source: /usr/local/cvs/OpenSees/SRC/api/elementAPI.h,v $

** Written: fmk
*/

#ifndef _eleAPI_G3
#define _eleAPI_G3

#include <string>
// #include <OPS_Globals.h>
#include <tcl.h>
// #include "TclModelBuilder.h"

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


// #ifdef __cplusplus
class TclSafeBuilder;
class UniaxialMaterial;
class NDMaterial;
class SectionForceDeformation;
class CrdTransf;
class FrictionModel;
class LimitCurve;
class Domain;
class FE_Datastore;

extern UniaxialMaterial* G3_getUniaxialMaterialInstance(Tcl_Interp*, int);
int G3_addUniaxialMaterial(Tcl_Interp*, UniaxialMaterial*);
extern Domain* G3_getDomain(Tcl_Interp*);
StaticAnalysis* G3_getStaticAnalysis(Tcl_Interp*);
int G3_setStaticAnalysis(Tcl_Interp*, StaticAnalysis*);
int G3_delStaticAnalysis(Tcl_Interp*);
DirectIntegrationAnalysis* G3_getTransientAnalysis(Tcl_Interp*);
int G3_setTransientAnalysis(Tcl_Interp*, DirectIntegrationAnalysis*);


StaticIntegrator* G3_getStaticIntegrator(Tcl_Interp*);
int G3_setStaticIntegrator(Tcl_Interp*, StaticIntegrator*);

TclSafeBuilder *G3_getSafeBuilder(Tcl_Interp *);

int G3_addTimeSeries(Tcl_Interp *, TimeSeries *);
TimeSeries* G3_getTimeSeries(Tcl_Interp *, int);

#endif // eleAPI_G3
