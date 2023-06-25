/* *****************************************************************************
  Copyright (c) 2015-2023, The Regents of the University of California (Regents).
  All rights reserved.

  Redistribution and use in source and binary forms, with or without 
  modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice, this
     list of conditions and the following disclaimer.
  2. Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer in the documentation
     and/or other materials provided with the distribution.

*************************************************************************** */
//
// Written: Minjie Zhu, cmp
//
#ifndef BasicAnalysisBulider_h
#define BasicAnalysisBulider_h

class Domain;
class G3_Table;
class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class EquiSolnAlgo;
class LinearSOE;
class EigenSOE;
class StaticIntegrator;
class TransientIntegrator;
class ConvergenceTest;
class StaticAnalysis;
class DirectIntegrationAnalysis;
class VariableTimeStepDirectIntegrationAnalysis;
class Integrator;

class BasicAnalysisBuilder
{
public:
    BasicAnalysisBuilder();
    BasicAnalysisBuilder(Domain* domain);
    ~BasicAnalysisBuilder();
/*
    int   tag_object(const char* type, int tag, void* obj);
    void* get_object(const char* type, int tag);
    void* pop_object(const char* type, int tag);
*/

    void set(ConstraintHandler* obj);
    void set(DOF_Numberer* obj);
    void set(EquiSolnAlgo* obj);
    void set(LinearSOE* obj);
    void set(Integrator* obj, int isstatic);
    void set(ConvergenceTest* obj);

    LinearSOE* getLinearSOE(int flag);
    
    Domain* getDomain(void);
    void newStaticAnalysis();
    int  newTransientAnalysis();
    int  setStaticAnalysis();
    int  setTransientAnalysis();
    void newEigenAnalysis(int typeSolver, double shift);
    int  getNumEigen() {return numEigen;};

    StaticAnalysis* getStaticAnalysis() {return theStaticAnalysis;}
    DirectIntegrationAnalysis* getTransientAnalysis() {return theTransientAnalysis;}
    VariableTimeStepDirectIntegrationAnalysis* getVariableTimeStepDirectIntegrationAnalysis() {
	return theVariableTimeStepTransientAnalysis;
    }

    EquiSolnAlgo*        getAlgorithm();
    StaticIntegrator*    getStaticIntegrator();
    TransientIntegrator* getTransientIntegrator();
    ConvergenceTest*     getConvergenceTest();

    int domainChanged(void);

    enum CurrentAnalysis {
      CURRENT_EMPTY_ANALYSIS,
      CURRENT_STATIC_ANALYSIS, 
      CURRENT_TRANSIENT_ANALYSIS
    } current_analysis;

    int analyze(int num_steps, double size_steps=0.0);
    int analyzeStatic(int num_steps);
    int analyzeTransient(int num_steps, double size_steps=0.0);

    void wipe();
    void resetStatic();
    void resetTransient();
    void resetAll();
    
private:
//  G3_Table* registry;
    enum CurrentAnalysis  CurrentAnalysisFlag = CURRENT_EMPTY_ANALYSIS;
    Domain                    *theDomain;
    ConstraintHandler 	      *theHandler;
    DOF_Numberer 	      *theNumberer;
    AnalysisModel 	      *theAnalysisModel;
    EquiSolnAlgo 	      *theAlgorithm;
    LinearSOE 		      *theSOE;
    EigenSOE 		      *theEigenSOE;
    StaticIntegrator          *theStaticIntegrator;
    TransientIntegrator       *theTransientIntegrator;
    ConvergenceTest           *theTest;
    StaticAnalysis            *theStaticAnalysis;
    DirectIntegrationAnalysis *theTransientAnalysis;
    VariableTimeStepDirectIntegrationAnalysis *theVariableTimeStepTransientAnalysis;

    int domainStamp;
    int numEigen;
};

#endif
