/* ****************************************************************** *\
|*    Opensees - Open System for Earthquake Engineering Simulation    *|
|*          Pacific Earthquake Engineering Research Center            *|
\* ****************************************************************** */

// $Date: 2021-10-30 $
// $Source: OpenSees/SRC/api/tclPackage/modelbuilder/TclSafeBuilder.h,v
// $

// Written: cmp
// Created: 10/21
//
// Description: This file contains the class definition for
// TclSafeBuilder. A TclSafeBuilder aims to be a threadsafe
// alternative to the TclBasicBuilder class. This class adds the commands to
// create the model for the standard models that can be generated using the
// elements released with the g3 framework.

#ifndef TCLSAFEBUILDER_H
#define TCLSAFEBUILDER_H

#include <ModelBuilder.h>
#include <MultiSupportPattern.h>

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

class TclSafeBuilder : public ModelBuilder {
public:
  TclSafeBuilder(Domain &theDomain, Tcl_Interp *interp, int ndm,
                         int ndf);
  ~TclSafeBuilder();

  int buildFE_Model(void);
  int getNDM(void) const;
  int getNDF(void) const;

  // methods needed for the truss and fiber-beam elements for
  // adding/getting uniaxial material objects
  // REMOVED    int addUniaxialMaterial(UniaxialMaterial &theMaterial);
  //            UniaxialMaterial *getUniaxialMaterial(int tag);

  // methods needed for the continuum elements and generic section
  // models to add/get ND material models

  //    int addNDMaterial(NDMaterial &theMaterial);
  //    NDMaterial *getNDMaterial(int tag);

  // methods needed for the nonlinear beam column elements to
  // add/get section objects
  /*
  int addSection(SectionForceDeformation &theSection);
  SectionForceDeformation *getSection(int tag);
  int addSectionRepres(SectionRepres &theSectionRepres);
  SectionRepres *getSectionRepres(int tag);
  */

  /*
  // methods needed for the yield surfaces
  int addYieldSurface_BC(YieldSurface_BC &theYS);
  YieldSurface_BC *getYieldSurface_BC(int tag);
  int addYS_EvolutionModel(YS_Evolution &theModel);
  YS_Evolution *getYS_EvolutionModel(int tag);
  int addPlasticMaterial(PlasticHardeningMaterial &theMaterial);
  PlasticHardeningMaterial *getPlasticMaterial(int tag);
  int addCyclicModel(CyclicModel &theModel); //!!
  CyclicModel *getCyclicModel(int tag); //!!
  int addDamageModel(DamageModel &theModel); //!!
  DamageModel *getDamageModel(int tag); //!!
  */

  // methods needed for the friction models
  // int addFrictionModel(FrictionModel &theFrnMdl);
  // FrictionModel *getFrictionModel(int tag);

  /* ----------------------------------------------- */
  Domain *getDomain();
  TclSafeBuilder *getBuilder();
  /* ----------------------------------------------- */

private:
  int ndm; // space dimension of the mesh
  int ndf; // number of degrees of freedom per node

  // TODO: change to std::map<>
  //    TaggedObjectStorage *theUniaxialMaterials;
  TaggedObjectStorage *theNDMaterials;
  TaggedObjectStorage *theSections;
  TaggedObjectStorage *theSectionRepresents;
  TaggedObjectStorage *theYieldSurface_BCs;
  TaggedObjectStorage *thePlasticMaterials;
  TaggedObjectStorage *theYS_EvolutionModels;
  TaggedObjectStorage *theCycModels; //!!
  //    TaggedObjectStorage *theDamageModels; //!!
  //    TaggedObjectStorage *theFrictionModels;
  TaggedObjectStorage *theLimitCurves; // MRL

  Domain *theTclDomain = 0;
  TclSafeBuilder *theTclBuilder = 0;
  int eleArgStart = 0;
  int nodeLoadTag = 0;
  int eleLoadTag = 0;

  // previously extern variables
  LoadPattern *theTclLoadPattern = 0;
  MultiSupportPattern *theTclMultiSupportPattern = 0;

protected:
  Tcl_Interp *theInterp;
};

#endif // TCLSAFEBUILDER_H

