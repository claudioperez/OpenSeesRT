/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** ****************************************************************** */
//
// Written: cmp
//
// Description: This file contains the implementaion of the
// TclCommand class.
//

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <tcl.h>
#include <G3_Logging.h>
#include <UniaxialMaterial.h>
#include <runtime/BasicModelBuilder.h>


//
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//

static Tcl_CmdProc TclCommand_useUniaxialMaterial;
static Tcl_CmdProc TclCommand_setStrainUniaxialMaterial;
static Tcl_CmdProc TclCommand_commitState;
static Tcl_CmdProc TclCommand_getStressUniaxialMaterial;
static Tcl_CmdProc TclCommand_getTangUniaxialMaterial;

const struct {const char*name; const Tcl_CmdProc*func;} command_table[] = {
  {"using",     TclCommand_useUniaxialMaterial       },
  {"with",      TclCommand_useUniaxialMaterial       },
  {"strain",    TclCommand_setStrainUniaxialMaterial },
  {"commit",    TclCommand_commitState               },
  {"stress",    TclCommand_getStressUniaxialMaterial },
  {"tangent",   TclCommand_getTangUniaxialMaterial   }
};


//
// THE FUNCTIONS INVOKED BY THE INTERPRETER
//
int
TclCommand_useUniaxialMaterial(ClientData clientData,
                                              Tcl_Interp *interp, int argc,
                                              TCL_Char **argv)
{
  // check number of arguments in command line
  if (argc < 4) {
    opserr << G3_ERROR_PROMPT << "bad arguments - want: using <obj-type> <obj-tag> {<operations>...}";
    return TCL_ERROR;
  }

  // get the tag form command line
  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "could not read obj-tag: using <obj-tag>?";
    return TCL_ERROR;
  }

  UniaxialMaterial *theMaterial = ((BasicModelBuilder*)clientData)->getUniaxialMaterial(argv[2]);
  // get the material from the modelbuilder with tag
  // and set the testing material to point to a copy of it
  if (theMaterial == nullptr) {
    opserr << G3_ERROR_PROMPT << "no material found with tag '" << argv[2] << "'.\n";
    return TCL_ERROR;

  } else {
    // theMaterial = theOrigMaterial->getCopy();
  }


  //
  // Add commands
  //
  const int ncmd = sizeof(command_table)/sizeof(command_table[0]);
  for (int i=0; i<ncmd; i++)
    Tcl_CreateCommand(interp,
                      command_table[i].name,
                      command_table[i].func,
                      (ClientData)theMaterial, 
                      nullptr);


  Tcl_Eval(interp, argv[3]);


  //
  //
  //
  Tcl_DeleteCommand(interp, "uniaxialTest");
  Tcl_DeleteCommand(interp, "strainUniaxialTest");
  Tcl_DeleteCommand(interp, "strain");
  Tcl_DeleteCommand(interp, "commit");
  Tcl_DeleteCommand(interp, "stressUniaxialTest");
  Tcl_DeleteCommand(interp, "tangUniaxialTest");

  return TCL_OK;
}

static int
TclCommand_setStrainUniaxialMaterial(ClientData clientData,
                                                    Tcl_Interp *interp,
                                                    int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  UniaxialMaterial* theMaterial = (UniaxialMaterial*)clientData;

  if (argc < 2) {
    opserr << G3_ERROR_PROMPT 
           << "bad arguments - want: strainUniaxialTest strain? <temp?>\n";
    return TCL_ERROR;
  }

  // get the tag from command line
  double strain;
  if (Tcl_GetDouble(interp, argv[1], &strain) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "could not read strain: strainUniaxialTest strain? "
              "<temp?>\n";
    return TCL_ERROR;
  }

  bool use_temp = false;
  bool commit = false;
  double temp = 0.0;
  if (argc > 2) {
    for (int i=2; i < argc; i++) {
      if (strcmp(argv[i], "-commit")==0){
        commit = true;
      } else if (Tcl_GetDouble(interp, argv[2], &temp) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "could not read strain: strainUniaxialTest strain? "
                  "<temp?>\n";
        return TCL_ERROR;
      }
    }
  }

  if (use_temp)
    theMaterial->setTrialStrain(strain, temp, 0.0);
  else
    theMaterial->setTrialStrain(strain);

  if (commit)
    theMaterial->commitState();

  return TCL_OK;
}

int TclCommand_commitState(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  UniaxialMaterial* theMaterial = (UniaxialMaterial*)clientData;
  theMaterial->commitState();
  return TCL_OK;
}

static int
TclCommand_getStressUniaxialMaterial(ClientData clientData,
                                                    Tcl_Interp *interp,
                                                    int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  UniaxialMaterial* theMaterial = (UniaxialMaterial*)clientData;
  double stress = 0.0;

  stress = theMaterial->getStress();
  char buffer[40];
  sprintf(buffer, "%.10e", stress);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);
  return TCL_OK;
}

static int
TclCommand_getTangUniaxialMaterial(ClientData clientData,
                                                  Tcl_Interp *interp, int argc,
                                                  TCL_Char **argv)
{
  assert(clientData != nullptr);
  UniaxialMaterial* theMaterial = (UniaxialMaterial*)clientData;

  double tangent = 0.0;

  tangent = theMaterial->getTangent();
  char buffer[40];
  sprintf(buffer, "%.10e", tangent);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);
  return TCL_OK;

}
