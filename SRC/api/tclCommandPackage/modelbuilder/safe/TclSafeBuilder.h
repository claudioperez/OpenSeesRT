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

#include <map>
#include <string>
// #include <ModelBuilder.h>
#include <TclBuilder.h>
#include <MultiSupportPattern.h>

class SectionForceDeformation;
class SectionRepres;
class NDMaterial;
class UniaxialMaterial;
class TaggedObjectStorage;
class YieldSurface_BC;
class YS_Evolution;
class PlasticHardeningMaterial;
class CyclicModel; //!!
class DamageModel;
class FrictionModel;
class TimeSeries;

#include <tcl.h>

class TclSafeBuilder : public TclBuilder {
public:
  TclSafeBuilder(Domain &theDomain, Tcl_Interp *interp, int ndm,
                         int ndf);
  ~TclSafeBuilder();

  using TclBuilder::buildFE_Model;
  // int buildFE_Model(void);

  int addTimeSeries(const std::string&, TimeSeries*);
  int addTimeSeries(TimeSeries*);
  TimeSeries* getTimeSeries(const std::string&);

  int incrNodalLoadTag(void);
  int decrNodalLoadTag(void);
  int getNodalLoadTag(void) const;

  LoadPattern* getEnclosingPattern(void) const;
  int setEnclosingPattern(LoadPattern*);

  // methods needed for the truss and fiber-beam elements for
  // adding/getting uniaxial material objects
  // REMOVED
  int addUniaxialMaterial(UniaxialMaterial *theMaterial);
  UniaxialMaterial *getUniaxialMaterial(int tag);

  // methods needed for the continuum elements and generic section
  // models to add/get ND material models

  //    int addNDMaterial(NDMaterial &theMaterial);
  NDMaterial *getNDMaterial(int tag);

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
  Domain *getDomain(void) const;
  TclSafeBuilder *getBuilder(void) const;
  /* ----------------------------------------------- */

private:
  int ndm; // space dimension of the mesh
  int ndf; // number of degrees of freedom per node

  // TODO: change to std::map<>
  TaggedObjectStorage *theUniaxialMaterials;
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
  LoadPattern *tclEnclosingPattern = 0;

  MultiSupportPattern *theTclMultiSupportPattern = 0;

  std::map<std::string, TimeSeries*> theTimeSeriesMap;

protected:
  Tcl_Interp *theInterp;
};

#endif // TCLSAFEBUILDER_H

