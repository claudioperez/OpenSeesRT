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
// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition for 
// WFFiberSection2d.h. WFFiberSection2d provides the abstraction of a 
// wide-flange section discretized by uniaxial fibers. 
//
#ifndef WFFiberSection2d_h
#define WFFiberSection2d_h

#include <FiberSection2d.h>
#include <Vector.h>
#include <Matrix.h>

#define SEC_TAG_WFFiberSection2d 1976

class NDMaterial;

class WFFiberSection2d : public FiberSection2d
{
 public:
  WFFiberSection2d(); 
  WFFiberSection2d(int tag, UniaxialMaterial &theMat,
	     double d, double tw, double bf, double tf,
	     int nfdw, int nftf);
  ~WFFiberSection2d();
  
  SectionForceDeformation *getCopy(void);
  
  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag = 0);

  const Vector& getStressResultantSensitivity(int gradIndex,
					      bool conditional);

 protected:
  
 private:
  double d;
  double tw;
  double bf;
  double tf;

  int nfdw;
  int nftf;
};

#endif
