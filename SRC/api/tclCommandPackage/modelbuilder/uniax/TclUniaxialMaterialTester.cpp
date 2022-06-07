/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** ****************************************************************** */

// Written: cmp
// Created: 03/01
//
// Description: This file contains the implementaion of the
// TclUniaxialMaterialTester class.
//

#include <stdlib.h>
#include <string.h>

#include <g3_api.h>
#include <ArrayOfTaggedObjects.h>
#include <UniaxialMaterial.h>
#include <TclUniaxialMaterialTester.h>


//
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//
typedef int (TclUniaxialTestCommand)(ClientData, Tcl_Interp*, int, TCL_Char**);

TclUniaxialTestCommand TclUniaxialMaterialTester_setUniaxialMaterial;
TclUniaxialTestCommand TclUniaxialMaterialTester_setStrainUniaxialMaterial;
TclUniaxialTestCommand TclUniaxialMaterialTester_commitState;
TclUniaxialTestCommand TclUniaxialMaterialTester_getStressUniaxialMaterial;
TclUniaxialTestCommand TclUniaxialMaterialTester_getTangUniaxialMaterial;

const struct {const char*name; const TclUniaxialTestCommand*func;} command_table[] = {
  {"using",     TclUniaxialMaterialTester_setUniaxialMaterial       },
  {"strain",    TclUniaxialMaterialTester_setStrainUniaxialMaterial },
  {"commit",    TclUniaxialMaterialTester_commitState               },
  {"stress",    TclUniaxialMaterialTester_getStressUniaxialMaterial },
  {"tangent",   TclUniaxialMaterialTester_getTangUniaxialMaterial   }
};

//
// CLASS CONSTRUCTOR & DESTRUCTOR
//

// constructor: the constructor will add certain commands to the interpreter
TclUniaxialMaterialTester::TclUniaxialMaterialTester(Domain &theDomain,
                                                     Tcl_Interp *interp,
                                                     int cTC)
    : TclSafeBuilder(theDomain, interp, 1, 1), theInterp(interp)
{
  const int ncmd = sizeof(command_table)/sizeof(command_table[0]);
  for (int i=0; i<ncmd; i++)
    Tcl_CreateCommand(interp,
                      command_table[i].name,
                      command_table[i].func,
                      (ClientData)NULL, 
                      NULL);
}

TclUniaxialMaterialTester::~TclUniaxialMaterialTester()
{


  Tcl_DeleteCommand(theInterp, "uniaxialTest");

  Tcl_DeleteCommand(theInterp, "strainUniaxialTest");
  Tcl_DeleteCommand(theInterp, "strain");
  Tcl_DeleteCommand(theInterp, "commit");

  Tcl_DeleteCommand(theInterp, "stressUniaxialTest");
  Tcl_DeleteCommand(theInterp, "tangUniaxialTest");
}

static UniaxialMaterial*
getUniaxialMaterial(Tcl_Interp *interp)
{
    return (UniaxialMaterial*)Tcl_GetAssocData(interp, "OPS::the_uniaxial_material", NULL);
}

//
// THE FUNCTIONS INVOKED BY THE INTERPRETER
//
int
TclUniaxialMaterialTester_setUniaxialMaterial(ClientData clientData,
                                              Tcl_Interp *interp, int argc,
                                              TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  UniaxialMaterial *theTestingUniaxialMaterial = getUniaxialMaterial(interp);

  // check number of arguments in command line
  if (argc < 4) {
    opserr << "WARNING bad arguments - want: using <obj-type> <obj-tag> {<operations>...}";
    return TCL_ERROR;
  }

  // get the matID form command line
  int matID;
  if (Tcl_GetInt(interp, argv[2], &matID) != TCL_OK) {
    opserr << "WARNING could not read obj-tag: using <obj-tag>?";
    return TCL_ERROR;
  }

  // delete the old testing material
  if (theTestingUniaxialMaterial != 0) {
    delete theTestingUniaxialMaterial;
    // theTestingUniaxialMaterial = 0;
    Tcl_SetAssocData(interp, "OPS::the_uniaxial_material", NULL, (ClientData)0);
  }

  // get the material from the modelbuilder with matID
  // and set the testing material to point to a copy of it
  UniaxialMaterial *theOrigMaterial = G3_getUniaxialMaterialInstance(rt, matID);
  if (theOrigMaterial == 0) {
    opserr << "WARNING no material found with tag '" << matID << "'.\n";
    return TCL_ERROR;

  } else {
    theTestingUniaxialMaterial = theOrigMaterial->getCopy();
    Tcl_SetAssocData(interp, "OPS::the_uniaxial_material", NULL, (ClientData)theTestingUniaxialMaterial);
  }

  if (argc > 3) Tcl_Eval(interp, argv[3]);

  return TCL_OK;
}

int
TclUniaxialMaterialTester_setStrainUniaxialMaterial(ClientData clientData,
                                                    Tcl_Interp *interp,
                                                    int argc, TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);

  if (argc < 2) {
    opserr
        << "WARNING bad arguments - want: strainUniaxialTest strain? <temp?>\n";
    return TCL_ERROR;
  }

  // get the matID form command line
  double strain;
  if (Tcl_GetDouble(interp, argv[1], &strain) != TCL_OK) {
    opserr << "WARNING could not read strain: strainUniaxialTest strain? "
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
        opserr << "WARNING could not read strain: strainUniaxialTest strain? "
                  "<temp?>\n";
        return TCL_ERROR;
      }
    }
  }

  // delete the old testing material
  UniaxialMaterial* theTestingUniaxialMaterial;
  if ((theTestingUniaxialMaterial=getUniaxialMaterial(interp))) {
    if (use_temp)
      theTestingUniaxialMaterial->setTrialStrain(strain, temp, 0.0);
    else
      theTestingUniaxialMaterial->setTrialStrain(strain);

    if (commit)
      theTestingUniaxialMaterial->commitState();
  }
  return TCL_OK;
}

int TclUniaxialMaterialTester_commitState(ClientData cd, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  UniaxialMaterial* theTestingUniaxialMaterial;
  if ((theTestingUniaxialMaterial=getUniaxialMaterial(interp))) {
      theTestingUniaxialMaterial->commitState();
  } else {
    opserr << "WARNING no active UniaxialMaterial - use uniaxialTest command\n";
    return TCL_ERROR;
  }
}

int
TclUniaxialMaterialTester_getStressUniaxialMaterial(ClientData clientData,
                                                    Tcl_Interp *interp,
                                                    int argc, TCL_Char **argv)
{
  double stress = 0.0;

  // delete the old testing material
  UniaxialMaterial* theTestingUniaxialMaterial;
  if ((theTestingUniaxialMaterial=getUniaxialMaterial(interp))) {
    stress = theTestingUniaxialMaterial->getStress();
    char buffer[40];
    sprintf(buffer, "%.10e", stress);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    //    sprintf(interp->result,"%.10e",stress);
    return TCL_OK;
  } else {
    opserr << "WARNING no active UniaxialMaterial - use uniaxialTest command\n";
    return TCL_ERROR;
  }
}

int
TclUniaxialMaterialTester_getTangUniaxialMaterial(ClientData clientData,
                                                  Tcl_Interp *interp, int argc,
                                                  TCL_Char **argv)
{
  double tangent = 0.0;
  // delete the old testing material
  UniaxialMaterial* theTestingUniaxialMaterial;
  if ((theTestingUniaxialMaterial=getUniaxialMaterial(interp))) {
    tangent = theTestingUniaxialMaterial->getTangent();
    char buffer[40];
    sprintf(buffer, "%.10e", tangent);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    //    sprintf(interp->result,"%.10e",tangent);
    return TCL_OK;
  } else {
    opserr << "WARNING no active UniaxialMaterial - use uniaxialTest command\n";
    return TCL_ERROR;
  }
}
