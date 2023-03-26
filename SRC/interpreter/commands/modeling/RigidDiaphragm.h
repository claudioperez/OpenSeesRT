/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Purpose: This file contains the class definition for RigidDiaphragm.
// RigidDiaphragm is a class which constructs MP_Constraint objects
// for a 3d Frame with a rigid diaphragm .. suitable for small
// displacement problems only.
//                                                                        
// $Revision: 1.3 $
// $Date: 2010-04-23 22:50:19 $
//
// File: ~/model/constraints/RigidDiaphragm.h
//
// Written: fmk 1/99

#ifndef RigidDiaphragm_h
#define RigidDiaphragm_h

class Domain;
class ID;

class RigidDiaphragm
{
  public:
    RigidDiaphragm(Domain &theDomain, int nodeR, ID &nodeC, 
		   int perpDirnToPlaneConstrained);
    virtual ~RigidDiaphragm(); 
  protected:   
  private:
};

#endif
