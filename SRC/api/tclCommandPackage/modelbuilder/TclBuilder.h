/* ****************************************************************** **
**    Opensees - Open System for Earthquake Engineering Simulation    **
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


// Description: This file contains the class definition for TclBuilder.
// A TclBuilder adds the commands to create the model for the standard
// models that can be generated using the elements released with the g3
// framework. currently these elements include:
//
// What: "@(#) TclBuilder.h, revA"

#ifndef TclBuilder_h
#define TclBuilder_h

#include <ModelBuilder.h>

class SectionForceDeformation;
class SectionRepres;
class NDMaterial;
class TaggedObjectStorage;
class YieldSurface_BC;
class YS_Evolution;
class PlasticHardeningMaterial;
class CyclicModel; //!!
class DamageModel;
class FrictionModel;

#include <tcl.h>
#include <g3_api.h>

class TclBuilder : public ModelBuilder {
public:
  TclBuilder(Domain &theDomain, int ndm, int ndf); // , Tcl_Interp *interp, int ndm, int ndf);
  ~TclBuilder();

  int buildFE_Model(void);
  int getNDM(void) const;
  int getNDF(void) const;

private:
  int ndm; // space dimension of the mesh
  int ndf; // number of degrees of freedom per node

protected:
  Tcl_Interp *theInterp;
};

#endif
