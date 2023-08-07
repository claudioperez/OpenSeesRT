/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Description: This file contains the class definition for
// BasicModelBuilder. A BasicModelBuilder aims to be a threadsafe
// alternative to the TclBasicBuilder class. This class adds the commands to
// create the model for the standard models that can be generated using the
// elements released with the g3 framework.
//
// Written: cmp
// Created: 10/21
//
#ifndef TCLSAFEBUILDER_H
#define TCLSAFEBUILDER_H

#include <tcl.h>
#include <string>
#include <unordered_map>
#include <TclBuilder.h>
#include <MultiSupportPattern.h>
#include "Storage/G3_TableIterator.h"

class SectionForceDeformation;
class SectionRepres;
class NDMaterial;
class UniaxialMaterial;
class TaggedObjectStorage;
class TimeSeries;
class G3_Runtime;
class CrdTransf;
class HystereticBackbone;
class G3_Table;


class BasicModelBuilder : public TclBuilder {
//
// CONSTRUCTORS / DESTRUCTORS
//
public:
  BasicModelBuilder(Domain &domain, Tcl_Interp *interp, int ndm, int ndf);
  ~BasicModelBuilder();
  using TclBuilder::buildFE_Model;
  typedef std::string key_t;
  template <typename ObjectType> class map_t
  : public std::unordered_map<key_t, ObjectType> {};

// Options
  void letClobber(bool option);
  bool canClobber();

  G3_TableIterator iterate(const char* partition);

//
// OBJECT CONTAINERS
// 
// Time series
private:
  G3_Table* registry = nullptr;
  map_t<TimeSeries*> m_TimeSeriesMap;
public:

  int   addRegistryObject(const char*, int tag, void* obj); 
  void* getRegistryObject(const char*, int tag); 

  int addTimeSeries(const std::string&, TimeSeries*);
  int addTimeSeries(TimeSeries*);
  TimeSeries* getTimeSeries(const key_t&);
  TimeSeries* getTimeSeries(int tag);
  LoadPattern* getEnclosingPattern(void);
  int setEnclosingPattern(LoadPattern*);
  int incrNodalLoadTag(void);
  int decrNodalLoadTag(void);
  int getNodalLoadTag(void);

  int addSP_Constraint(int axisDirn, 
         double axisValue, 
         const ID &fixityCodes, 
         double tol=1e-10);

// Coordinate Transformations
private: map_t<CrdTransf*> m_CrdTransfMap;
public:  int addCrdTransf(CrdTransf *);
         int addCrdTransf(const key_t, CrdTransf*);
         CrdTransf *getCrdTransf(int tag);
         CrdTransf *getCrdTransf(const key_t&);

// Uniaxial materials
private: map_t<UniaxialMaterial*> m_UniaxialMaterialMap;
public:  int  addUniaxialMaterial(UniaxialMaterial &theMaterial);
         int  addUniaxialMaterial(const std::string&, UniaxialMaterial &);
         int  addUniaxialMaterial(UniaxialMaterial *theMaterial);
         UniaxialMaterial *getUniaxialMaterial(int tag);
         UniaxialMaterial *getUniaxialMaterial(const std::string &);

// Backbone materials
private: map_t<HystereticBackbone*> m_HystereticBackboneMap;
public:  int addHystereticBackbone(HystereticBackbone &theMaterial);
         int addHystereticBackbone(const std::string&, HystereticBackbone &);
         // HystereticBackbone *getHystereticBackbone(int tag);
         HystereticBackbone *getHystereticBackbone(const std::string &);

// Multi-dimensional materials
private: map_t<NDMaterial*> m_NDMaterialMap;
public:  int addNDMaterial(NDMaterial &theMaterial);
         int addNDMaterial(const std::string&, NDMaterial &);
         NDMaterial *getNDMaterial(int tag);
         NDMaterial *getNDMaterial(const std::string &);

// Cross sections
private: map_t<SectionForceDeformation*> m_SectionForceDeformationMap;
         map_t<SectionRepres*          > m_SectionRepresMap;
public:  int addSection(SectionForceDeformation &theSection);
         int addSection(const std::string&, SectionForceDeformation &);
         SectionForceDeformation *getSection(int tag);
         SectionForceDeformation *getSection(const std::string &);
         int addSectionRepres(SectionRepres &theSectionRepres);
         int addSectionRepres(const std::string &, SectionRepres &);
         SectionRepres *getSectionRepres(int tag);
         SectionRepres *getSectionRepres(const std::string&);

//
// OTHER METHODS
//
  Domain *getDomain(void) const;
  BasicModelBuilder *getBuilder(void) const;

private:
  int ndm; // space dimension of the mesh
  int ndf; // number of degrees of freedom per node

  G3_Runtime *m_runtime = nullptr;
  Domain *theTclDomain = 0;
  BasicModelBuilder *theTclBuilder = 0;
  int eleArgStart = 0;
  int nodeLoadTag = 0;
  int eleLoadTag = 0;

  // Options
  bool no_clobber = true;

// previously extern variables
  LoadPattern *tclEnclosingPattern = nullptr;
  MultiSupportPattern *theTclMultiSupportPattern = nullptr;

protected:
  Tcl_Interp *theInterp;

};

#endif // TCLSAFEBUILDER_H

