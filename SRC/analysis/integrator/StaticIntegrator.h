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
// Description: This file contains the class definition for StaticIntegrator.
// StaticIntegrator is an algorithmic class for setting up the finite element
// equations for a static analysis and for Incrementing the nodal displacements
// with the values in the soln vector to the LinearSOE object. 
//
// Written: fmk 
// Created: 11/96
// File: ~/analysis/integrator/StaticIntegrator.h
//
#ifndef StaticIntegrator_h
#define StaticIntegrator_h

#include <IncrementalIntegrator.h>

class FE_Element;
class DOF_Group;

enum class ReferencePattern {
  Point,
  Body,
  Full
};



class StaticIntegrator : public IncrementalIntegrator
{
  public:
    StaticIntegrator(int classTag);    

    virtual ~StaticIntegrator();

    virtual int newStep() = 0;

    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    virtual int formUnbalance() final;
    virtual int formEleTangent(FE_Element *theEle);
    virtual int formEleResidual(FE_Element *theEle)   final;
    virtual int formNodTangent(DOF_Group *theDof)     final;
    virtual int formNodUnbalance(DOF_Group *theDof)   final;    
    virtual int formEleTangentSensitivity(FE_Element *theEle,int gradNumber); 

  protected:
    enum class ResidualType {
      StaticUnbalance,
      StaticSensitivity,
      TransientUnbalance,
//    TransientSensitivity
    };
    void setResidualType(ResidualType tag) {residualType = tag;}
    ResidualType  getResidualType() { return residualType;}

    void setGradIndex(int index) {gradIndex = index;}
    int  getGradIndex() { return gradIndex;}

 
  private:
    ResidualType residualType;
    int gradIndex;
};

#endif

