//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for BasicModelBuilder.
// A BasicModelBuilder adds the commands to create the model for the standard
// models that can be generated using the elements released with the g3
// framework.
//
// TODO:
// - Remove all *Map.at
//   - return from registry
//   - handle find failures consistently
//
// Written: cmp
//
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <initializer_list>

#include <modeling/commands.h>

#include <runtimeAPI.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Domain.h>

#include <CrdTransf.h>
#include <SectionForceDeformation.h>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <runtime/BasicModelBuilder.h>
#include <MultiSupportPattern.h>

#include <TimeSeries.h>

#include <tcl.h> // For TCL_OK/ERROR

//
// CLASS CONSTRUCTOR & DESTRUCTOR
//
BasicModelBuilder::BasicModelBuilder(Domain &domain, Tcl_Interp *interp, 
                                     int NDM, int NDF)
    : ndm(NDM), ndf(NDF), theInterp(interp),
      section_builder_is_set(false),
      theDomain(&domain),
      tclEnclosingPattern(nullptr),
      next_node_load(0),
      next_elem_load(0)

{
  static int ncmd = sizeof(tcl_char_cmds)/sizeof(char_cmd);

  Tcl_CreateCommand(interp, "wipe", TclCommand_wipeModel, (ClientData)this, nullptr);

  for (int i = 0; i < ncmd; i++)
    Tcl_CreateCommand(interp, 
        tcl_char_cmds[i].name, 
        tcl_char_cmds[i].func, 
        (ClientData) this, nullptr);
 
  tclEnclosingPattern = nullptr;
  // theTclMultiSupportPattern = 0;

  Tcl_SetAssocData(interp, "OPS::theTclBuilder", NULL, (ClientData)this);
  Tcl_SetAssocData(interp, "OPS::theBasicModelBuilder", NULL, (ClientData)this);
  Tcl_SetAssocData(interp, "OPS::theTclDomain", NULL, (ClientData)&theDomain);

}

BasicModelBuilder::~BasicModelBuilder()
{

  for (auto& [part, val] : m_registry ) {
    for (auto& [tag, obj] : val)
      delete obj;
  }

  // set the pointers to 0
  theDomain = nullptr;
//theTclBuilder = nullptr;
  tclEnclosingPattern = nullptr;

  static int ncmd = sizeof(tcl_char_cmds)/sizeof(char_cmd);
  for (int i = 0; i < ncmd; i++)
    Tcl_DeleteCommand(theInterp, tcl_char_cmds[i].name);
}


int
BasicModelBuilder::buildFE_Model() {return 0;}

int
BasicModelBuilder::getNDM() const {return ndm;}

int
BasicModelBuilder::getNDF() const {return ndf;}

LoadPattern*
BasicModelBuilder::getCurrentLoadPattern() 
{
  return m_current_load_pattern;
}


void
BasicModelBuilder::letClobber(bool let_clobber) {
  no_clobber = !let_clobber;
}

bool
BasicModelBuilder::canClobber() {
  return !no_clobber;
}

int BasicModelBuilder::incrNodalLoadTag(){return ++next_node_load;};
int BasicModelBuilder::decrNodalLoadTag(){return --next_node_load;};
int BasicModelBuilder::getNodalLoadTag() {return   next_node_load;};

int
BasicModelBuilder::addSP_Constraint(int axisDirn, double axisValue, const ID &fixityCodes, double tol)
{
  return theDomain->addSP_Constraint(axisDirn, axisValue, fixityCodes, tol);
}

LoadPattern *
BasicModelBuilder::getEnclosingPattern()
{
  return tclEnclosingPattern;
}

int
BasicModelBuilder::setEnclosingPattern(LoadPattern* pat)
{
  tclEnclosingPattern = pat;
  return 1;
}

int
BasicModelBuilder::getCurrentSectionBuilder(int& tag)
{
  if (section_builder_is_set) {
    tag = current_section_builder;
    return  0;
  } else
    return -1;
}

void 
BasicModelBuilder::setCurrentSectionBuilder(int tag)
{
  section_builder_is_set   = true;
  current_section_builder  = tag;
}

Domain *
BasicModelBuilder::getDomain() const 
{
  return theDomain;
}

#if 0
BasicModelBuilder *
BasicModelBuilder::getBuilder() const {
  return theTclBuilder;
}
#endif

int 
BasicModelBuilder::printRegistry(const char *partition, OPS_Stream& stream, int flag) const 
{
    int count = 0;
    auto iter = m_registry.find(partition);
    if (iter == m_registry.end()) {
      return count;
    }

    for (auto const& [key, val] : iter->second) {
      if (count != 0)
        stream << ",\n";

      val->Print(stream, flag);
      count++;
    }

    return count;
}

void* 
BasicModelBuilder::getRegistryObject(const char* partition, int tag) const
{

  auto iter = m_registry.find(std::string{partition});
  if (iter == m_registry.end()) {
    opserr << "No objects of type \"" << partition << "\" have been created.\n";
    return nullptr;
  }

  auto iter_objs = iter->second.find(tag) ;
  if (iter_objs == iter->second.end()) {
    opserr << "No object with tag \"" << tag << "\"in partition \"" << partition << "\"\n";
    return nullptr;
  }

  return (void*)iter_objs->second;

}

int
BasicModelBuilder::addRegistryObject(const char* partition, int tag, void *obj)
{
  // TODO: Change void* obj to TaggedObject*
  m_registry[std::string{partition}][tag] = (TaggedObject*)obj;
  return TCL_OK;
}

