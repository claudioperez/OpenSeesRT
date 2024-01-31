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
// Fiber. Fiber is an abstract base class and thus no objects of 
// it's type can be instatiated. It has pure virtual functions which
// must be implemented in it's derived classes.
//
// Written: Remo Magalhaes de Souza
// Created: 10/98
//
#ifndef Fiber_h
#define Fiber_h
//
#include <TaggedObject.h>
#include <MovableObject.h>

class Matrix;
class ID;
class UniaxialMaterial;
class NDMaterial;
class Information;
class Response;

class Fiber : public TaggedObject, public MovableObject
{
  public:
    Fiber (int tag, int classTag, double y, double z, double area);
    virtual ~Fiber();

    virtual int    setTrialFiberStrain(const Vector &vs)=0;
    virtual Vector &getFiberStressResultants(void) =0;
    virtual Matrix &getFiberTangentStiffContr(void) =0;
	
    virtual int    commitState(void)=0;
    virtual int    revertToLastCommit(void)=0;    
    virtual int    revertToStart(void)=0;
    
    virtual Fiber *getCopy(void) = 0;
    virtual int getOrder(void) = 0;
    virtual const ID &getType(void) = 0;

    virtual Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    virtual int getResponse(int responseID, Information &info);

    double getArea(void) const {
      return area;
    };
    void getFiberLocation(double &yLoc, double &zLoc) const {
      yLoc = loc_y; //  yLoc = -y;
      zLoc = loc_z; //  zLoc = 0.0;
    }

    virtual double getd(void) =0;

    virtual UniaxialMaterial *getMaterial(void) {return 0;}
    virtual NDMaterial       *getNDMaterial(void) {return 0;}

    virtual const Vector &getFiberSensitivity(int gradNumber, bool cond);
    virtual int commitSensitivity(const Vector &dedh, int gradNumber,
				  int numGrads);


 protected:
    Vector *sDefault;
    Matrix *fDefault;
    double area;
    const double loc_y, 
                 loc_z = 0.0;

 private:
};

#endif
