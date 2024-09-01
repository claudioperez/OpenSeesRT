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
// Purpose: This file contains the class definition for Load.
// Load is an abstract class. A Load object is used to add load
// to the model. 
//
// File: ~/domain/load/Load.h
//
// Written: fmk 
// Created: 11/96
//
#ifndef Load_h
#define Load_h

#include <ID.h>
#include <TaggedObject.h>
#include <MovableObject.h>
class Domain;

class Load : public TaggedObject, public MovableObject
{
  public:
    Load(int tag, int classTag);

    virtual ~Load();

    // pure virtual functions
    virtual void applyLoad(double loadfactor) =0;
    virtual void applyLoadSensitivity(double loadfactor) {return;}
    
    virtual void setLoadPatternTag(int loadPaternTag);
    virtual int  getLoadPatternTag(void) const;

#if 1
    virtual void setDomain(Domain *theDomain);
    virtual Domain *getDomain(void) const final;
  private:
    Domain* theDomain;
#endif

  protected:
	
  private:
    int loadPatternTag;
};

#endif
