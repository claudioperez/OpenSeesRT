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

#include <typeinfo>
// #include <tcl.h>
#include <string>
#include <unordered_map>
#include <runtime/modelbuilder/TclBuilder.h>


#include <TaggedObject.h>
class MultiSupportPattern;
class G3_Runtime;
class OPS_Stream;
class ID;
struct Tcl_Interp;

class BasicModelBuilder : public TclBuilder {
public:
//
// CONSTRUCTORS / DESTRUCTORS
//
  BasicModelBuilder(Domain &domain, Tcl_Interp *interp, int ndm, int ndf);
  ~BasicModelBuilder();

  using TclBuilder::buildFE_Model;


// Options
  void letClobber(bool option);
  bool canClobber();


  int   addRegistryObject(const char*, int tag, void* obj); 

  template<class T> int addTypedObject(int tag, T* obj) {
    return addRegistryObject(typeid(T).name(), tag, obj);
  }

  template<class T> int addTaggedObject(T& obj) {
    int tag = obj.getTag();
    m_registry[typeid(T).name()][tag] = &obj;
    return addRegistryObject(typeid(T).name(), tag, &obj);
  }

//void printRegistry(const char* partition, OPS_Stream& stream, int flag) const;
  template <class T>
  void printRegistry(OPS_Stream& stream, int flag) const 
  {
    auto partition = typeid(T).name();

    printRegistry(partition, stream, flag);

  }


  void* getRegistryObject(const char*, int tag) const;

  template<class T> T* getTypedObject(int tag) {
    return (T*)getRegistryObject(typeid(T).name(), tag);
  }

  LoadPattern* getEnclosingPattern(void);
  int setEnclosingPattern(LoadPattern*);
  int incrNodalLoadTag(void);
  int decrNodalLoadTag(void);
  int getNodalLoadTag(void);

  int addSP_Constraint(int axisDirn, 
         double axisValue, 
         const ID &fixityCodes, 
         double tol=1e-10);

//
// OTHER METHODS
//
  Domain *getDomain(void) const;
  BasicModelBuilder *getBuilder(void) const;

protected:
  Tcl_Interp *theInterp;

// 
private:
  void printRegistry(const char *, OPS_Stream& stream, int flag) const ;


  int ndm; // space dimension of the mesh
  int ndf; // number of degrees of freedom per node

  G3_Runtime *m_runtime = nullptr;
  Domain *theTclDomain = 0;
  BasicModelBuilder *theTclBuilder = nullptr;
  int eleArgStart = 0;
  int nodeLoadTag = 0;
  int eleLoadTag = 0;

  // Options
  bool no_clobber = true;

// previously extern variables
  LoadPattern *tclEnclosingPattern = nullptr;
  MultiSupportPattern *theTclMultiSupportPattern = nullptr;

// OBJECT CONTAINERS
  std::unordered_map<std::string, std::unordered_map<int, TaggedObject*>> m_registry;

};

#endif

