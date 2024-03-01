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
// Description: This file contains the implementation for the
// Fiber class. Fiber provides the abstraction of a section fiber.
//
// Written: Remo Magalhaes de Souza
// Created: 10/98
// Revision: 
//
#include <Fiber.h>
#include <Matrix.h>
#include <Vector.h>

// constructor:
Fiber::Fiber(int tag, int classTag, double y, double z, double area):
  TaggedObject(tag), MovableObject(classTag),
  loc_y(y), loc_z(z), area(area),
  sDefault(0), fDefault(0)
{
}

// destructor:
Fiber::~Fiber()
{
  if (sDefault != 0)
    delete sDefault;
  if (fDefault != 0)
    delete fDefault;
}

Response*
Fiber::setResponse(const char **argv, int argc, OPS_Stream &s)
{
  return 0;
}

int
Fiber::getResponse(int responseID, Information &info)
{
  return -1;
}

const Vector&
Fiber::getFiberSensitivity(int gradNumber, bool cond)
{
  if (sDefault == 0)
    sDefault = new Vector (this->getOrder());
  return *sDefault;

}

int
Fiber::commitSensitivity(const Vector &dedh, int gradNumber,
			 int numGrads)
{
  return -1;
}
