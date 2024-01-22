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
//
// Description: This file contains the class definition for 
// DirectIntegrationAnalysis. DirectIntegrationAnalysis is a 
// subclass of TransientAnalysis. It is used to perform a 
// dynamic analysis on the FE\_Model using a direct integration scheme.
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
#ifndef DirectIntegrationAnalysis_h
#define DirectIntegrationAnalysis_h
//
// #include <Analysis.h>

class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class TransientIntegrator;
class LinearSOE;
class EquiSolnAlgo;
class ConvergenceTest;
class EigenSOE;
class Domain;

// class DirectIntegrationAnalysis: public TransientAnalysis
class DirectIntegrationAnalysis // : public Analysis
{
  public:
    DirectIntegrationAnalysis(Domain &theDomain, 
			      ConstraintHandler &theHandler,
			      DOF_Numberer &theNumberer,
			      AnalysisModel &theModel,
			      EquiSolnAlgo &theSolnAlgo,		   
			      LinearSOE &theSOE,
			      TransientIntegrator &theIntegrator,
			      ConvergenceTest *theTest = 0,
			      int numSubLevels = 0,
			      int numSubSteps = 0);


    virtual ~DirectIntegrationAnalysis();

    void clearAll(void);	    
    
    int analyze(int numSteps, double dT);
    int analyzeStep(double dT);
    int analyzeSubLevel(int level, double dT);
    int eigen(int numMode, bool generlzed = true, bool findSmallest = true);
    int initialize(void);
    int domainChanged(void);

    int setNumberer(DOF_Numberer &theNumberer);    
#if 0
    int setAlgorithm(EquiSolnAlgo &theAlgorithm);
    int setIntegrator(TransientIntegrator &theIntegrator);
#endif
    int setLinearSOE(LinearSOE &theSOE); 
    int setConvergenceTest(ConvergenceTest &theTest);
    int setEigenSOE(EigenSOE &theSOE);
    
    int checkDomainChange(void);

    EquiSolnAlgo        *getAlgorithm(void);
    TransientIntegrator *getIntegrator(void);
    ConvergenceTest     *getConvergenceTest(void); 
    //AnalysisModel       *getModel(void) ;

  protected: 
    Domain              *theDomain;
    AnalysisModel 	*theAnalysisModel;

  private:
    ConstraintHandler 	*theConstraintHandler;    
    DOF_Numberer 	*theDOF_Numberer;
    EquiSolnAlgo 	*theAlgorithm;
    LinearSOE 		*theSOE;
    EigenSOE 		*theEigenSOE;
    TransientIntegrator *theIntegrator;
    ConvergenceTest     *theTest;

    int domainStamp;
    int numSubLevels;
    int numSubSteps;


};

#endif

