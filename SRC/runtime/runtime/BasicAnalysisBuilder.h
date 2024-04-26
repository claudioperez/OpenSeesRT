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
// BasicAnalysisBuilder is an aggregate class which manages the analysis objects:
//
// - LinearSOE
// - Domain                    *theDomain;
// - ConstraintHandler 	      *theHandler;
// - DOF_Numberer 	      *theNumberer;
// - AnalysisModel 	      *theAnalysisModel;
// - EquiSolnAlgo 	      *theAlgorithm;
// - EigenSOE 		      *theEigenSOE;
// - StaticIntegrator          *theStaticIntegrator;
// - TransientIntegrator       *theTransientIntegrator;
// - ConvergenceTest           *theTest;
//
// The BasicAnalysisBuilder assumes responsibility for
// deleting these objects, but ownership of the SOE may be
// given up.
//
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
class VariableTimeStepDirectIntegrationAnalysis;
class Integrator;

class BasicAnalysisBuilder
{
public:
    BasicAnalysisBuilder(Domain* domain);
    ~BasicAnalysisBuilder();

//  enum NoDelete {
//    StaticIntegrator    = 1<<0,
//    TransientIntegrator = 1<<1,
//    LinearSOE           = 1<<2,
//  };

    enum CurrentAnalysis {
      EMPTY_ANALYSIS,
      STATIC_ANALYSIS, 
      TRANSIENT_ANALYSIS
    };

    void set(ConstraintHandler* obj);
    void set(DOF_Numberer* obj);
    void set(EquiSolnAlgo* obj);
    void set(LinearSOE*  obj, bool free=true);
    void set(StaticIntegrator& obj);
    void set(TransientIntegrator& obj, bool free=true);
    void set(ConvergenceTest* obj);
    void set(EigenSOE& obj);

    LinearSOE* getLinearSOE();

    Domain* getDomain(void);
    int initialize(void);

    int  newTransientAnalysis();
    int  setStaticAnalysis();
    int  setTransientAnalysis();

    //   Eigen
    void newEigenAnalysis(int typeSolver, double shift);
    int  eigen(int numMode, bool generalized, bool findSmallest);
    int  getNumEigen() {return numEigen;};

    void formUnbalance();

    VariableTimeStepDirectIntegrationAnalysis* getVariableTimeStepDirectIntegrationAnalysis() {
	return theVariableTimeStepTransientAnalysis;
    }

    EquiSolnAlgo*        getAlgorithm();
    StaticIntegrator*    getStaticIntegrator();
    TransientIntegrator* getTransientIntegrator();
    ConvergenceTest*     getConvergenceTest();

    int domainChanged(void);

    int analyze(int num_steps, double size_steps=0.0);
    int analyzeStatic(int num_steps);
    
    int analyzeTransient(int numSteps, double dT);
    int analyzeStep(double dT);
    int analyzeSubLevel(int level, double dT);

    void wipe();

    
    enum CurrentAnalysis  CurrentAnalysisFlag = EMPTY_ANALYSIS;

private:
    void setLinks(CurrentAnalysis flag = EMPTY_ANALYSIS);
    void fillDefaults(enum CurrentAnalysis flag);

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
    VariableTimeStepDirectIntegrationAnalysis *theVariableTimeStepTransientAnalysis;

    int domainStamp;
    int numEigen = 0;

    int numSubLevels = 0;
    int numSubSteps  = 0;

    bool freeSOE = true;
    bool freeTI  = true;

};

#endif
