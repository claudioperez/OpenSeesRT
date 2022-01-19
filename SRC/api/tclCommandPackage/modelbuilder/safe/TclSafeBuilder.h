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
class LimitCurve;
class DamageModel;
class FrictionModel;
class TimeSeries;
class G3_Runtime;

#include <tcl.h>

class TclSafeBuilder : public TclBuilder {
//
// CONSTRUCTORS / DESTRUCTORS
//
public:
  TclSafeBuilder(Domain &domain, Tcl_Interp *interp, int ndm, int ndf);
  ~TclSafeBuilder();
  using TclBuilder::buildFE_Model;


//
// OBJECT CONTAINERS
// 
   // Time series
private:
  std::map<std::string, TimeSeries*> m_TimeSeriesMap;
public:
  int addTimeSeries(const std::string&, TimeSeries*);
  int addTimeSeries(TimeSeries*);
  TimeSeries* getTimeSeries(const std::string&);
  LoadPattern* getEnclosingPattern(void) const;
  int setEnclosingPattern(LoadPattern*);
  int incrNodalLoadTag(void);
  int decrNodalLoadTag(void);
  int getNodalLoadTag(void) const;

  // Uniaxial materials
private: std::map<std::string, UniaxialMaterial*> m_UniaxialMaterialMap;
public:  int  addUniaxialMaterial(UniaxialMaterial &theMaterial);
         int  addUniaxialMaterial(const std::string&, UniaxialMaterial &);
         int  addUniaxialMaterial(UniaxialMaterial *theMaterial);
         UniaxialMaterial *getUniaxialMaterial(int tag);
         UniaxialMaterial *getUniaxialMaterial(const std::string &);


  // Multi-dimensional materials
private: std::map<std::string,  NDMaterial*> m_NDMaterialMap;
public:  int         addNDMaterial(NDMaterial &theMaterial);
         int         addNDMaterial(const std::string&, NDMaterial &);
         NDMaterial *getNDMaterial(int tag);
         NDMaterial *getNDMaterial(const std::string &);

  // Cross sections
private: std::map<std::string, SectionForceDeformation*> m_SectionForceDeformationMap;
         std::map<std::string, SectionRepres*          > m_SectionRepresMap;
public:  int addSection(SectionForceDeformation &theSection);
         int addSection(const std::string&, SectionForceDeformation &);
         SectionForceDeformation *getSection(int tag);
         SectionForceDeformation *getSection(const std::string &);
         int addSectionRepres(SectionRepres &theSectionRepres);
         int addSectionRepres(const std::string &, SectionRepres &);
         SectionRepres *getSectionRepres(int tag);
         SectionRepres *getSectionRepres(const std::string&);

  // methods needed for the yield surfaces
/*
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

  // methods needed for the friction models
  int addFrictionModel(FrictionModel &theFrnMdl);
  FrictionModel *getFrictionModel(int tag);
*/

  /* ----------------------------------------------- */
  Domain *getDomain(void) const;
  TclSafeBuilder *getBuilder(void) const;
  /* ----------------------------------------------- */

private:
  int ndm; // space dimension of the mesh
  int ndf; // number of degrees of freedom per node

  // TODO: change to std::map<>
  TaggedObjectStorage *theNDMaterials;
  // TaggedObjectStorage *theSections;
  TaggedObjectStorage *theSectionRepresents;
  TaggedObjectStorage *theYieldSurface_BCs;
  TaggedObjectStorage *thePlasticMaterials;
  TaggedObjectStorage *theYS_EvolutionModels;
  TaggedObjectStorage *theCycModels; //!!
  //    TaggedObjectStorage *theDamageModels; //!!
  //    TaggedObjectStorage *theFrictionModels;
  std::map<std::string, LimitCurve*> theLimitCurves; // MRL


  Domain *theTclDomain = 0;
  G3_Runtime *m_runtime = nullptr;
  TclSafeBuilder *theTclBuilder = 0;
  int eleArgStart = 0;
  int nodeLoadTag = 0;
  int eleLoadTag = 0;

  // previously extern variables
  LoadPattern *tclEnclosingPattern = 0;

  MultiSupportPattern *theTclMultiSupportPattern = 0;

protected:
  Tcl_Interp *theInterp;
};

#endif // TCLSAFEBUILDER_H

