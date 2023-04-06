/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Written: fmk 12/99
//
// Purpose: This file contains the class definition for RigidRod.
// RigidRod is a class which constructs MP_Constraint objects
// for a rigid rod, all translational dof are constrained to be equal
// at the retained and constarined nodes.

#ifndef RigidRod_h
#define RigidRod_h

class Domain;
class ID;

class RigidRod
{
  public:
    RigidRod(Domain &theDomain, int nodeR, int nodeC);
    virtual ~RigidRod();
    
  protected:
    
  private:
};

#endif
