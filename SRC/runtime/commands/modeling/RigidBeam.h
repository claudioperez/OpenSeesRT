/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Purpose: This file contains the class definition for RigidBeam.
// RigidBeam is a class which constructs an MP_Constraint object
// between two nodes which is similar to rigid beam
//
// Written: fmk 12/99
//                                                                        
// File: ~/model/constraints/RigidBeam.h
//

#ifndef RigidBeam_h
#define RigidBeam_h

class Domain;
class ID;

class RigidBeam
{
  public:
    RigidBeam(Domain &theDomain, int nodeR, int nodeC);
    virtual ~RigidBeam();
    
  protected:
    
  private:
};

#endif
