//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#include <tcl.h>
#include <string>
#include <unordered_map>

Tcl_CmdProc TclCommand_newElasticMaterial;

static
std::unordered_map<std::string, Tcl_CmdProc*> MaterialLibrary = {
  {"ElasticIsotropic",          TclCommand_newElasticMaterial},
};
