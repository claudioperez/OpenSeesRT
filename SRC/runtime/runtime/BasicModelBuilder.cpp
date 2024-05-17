/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
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
BasicModelBuilder::BasicModelBuilder(Domain &theDomain, Tcl_Interp *interp, int NDM,
                               int NDF)
    : TclBuilder(theDomain, NDM, NDF), theInterp(interp)
{
  static int ncmd = sizeof(tcl_char_cmds)/sizeof(char_cmd);

  Tcl_CreateCommand(interp, "wipe", TclCommand_wipeModel, (ClientData)this, nullptr);

  for (int i = 0; i < ncmd; i++)
    Tcl_CreateCommand(interp, 
        tcl_char_cmds[i].name, 
        tcl_char_cmds[i].func, 
        (ClientData) this, nullptr);
 
  theTclBuilder = this;
  theTclDomain = &theDomain;
  tclEnclosingPattern = nullptr;
  // theTclMultiSupportPattern = 0;

  nodeLoadTag = 0;
  eleArgStart = 0;

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
  theTclDomain = nullptr;
  theTclBuilder = nullptr;
  tclEnclosingPattern = nullptr;

  static int ncmd = sizeof(tcl_char_cmds)/sizeof(char_cmd);
  for (int i = 0; i < ncmd; i++)
    Tcl_DeleteCommand(theInterp, tcl_char_cmds[i].name);
}

//
// CLASS METHODS
//
void
BasicModelBuilder::letClobber(bool let_clobber) {
  no_clobber = !let_clobber;
}

bool
BasicModelBuilder::canClobber() {
  return !no_clobber;
}

int BasicModelBuilder::incrNodalLoadTag(void){return ++nodeLoadTag;};
int BasicModelBuilder::decrNodalLoadTag(void){return --nodeLoadTag;};
int BasicModelBuilder::getNodalLoadTag(void) {return   nodeLoadTag;};

int
BasicModelBuilder::addSP_Constraint(int axisDirn, double axisValue, const ID &fixityCodes, double tol)
{
  return theTclDomain->addSP_Constraint(axisDirn, axisValue, fixityCodes, tol);
}

LoadPattern *
BasicModelBuilder::getEnclosingPattern(void)
{
  return tclEnclosingPattern;
}

int
BasicModelBuilder::setEnclosingPattern(LoadPattern* pat)
{
  tclEnclosingPattern = pat;
  return 1;
}

Domain *
BasicModelBuilder::getDomain(void) const 
{
  return theTclDomain;
}

BasicModelBuilder *
BasicModelBuilder::getBuilder(void) const {
  return theTclBuilder;
}

void 
BasicModelBuilder::printRegistry(const char *partition, OPS_Stream& stream, int flag) const 
{
    auto iter = m_registry.find(partition);
    if (iter == m_registry.end()) {
      // opserr << "No objects of type \"" << partition << "\" have been created.\n";
      return;// nullptr;
    }

    bool first = true;
    for (auto const& [key, val] : iter->second) {
      if (!first)
        stream << ",\n";

      val->Print(stream, flag);

      first = false;
    }
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

