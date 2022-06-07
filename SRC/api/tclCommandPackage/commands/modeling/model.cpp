/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Description: This file contains the function myCommands().
// myCommands() is called in g3AppInit() - all new user commands
// are to be placed in here.
//
#include <g3_api.h>
#include <Domain.h>
// #include "TclBasicBuilder.h"
#include "TclUniaxialMaterialTester.h"
#include "TclPlaneStressMaterialTester.h"
#include "modelbuilder/safe/TclSafeBuilder.h"
#include "modelbuilder/sect/TclSectionTestBuilder.h"

#include <g3_api.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern ModelBuilder *theBuilder;

#ifdef _PARALLEL_PROCESSING
#include <PartitionedDomain.h>
extern PartitionedDomain theDomain;
#else
extern Domain theDomain;
#endif

int specifyModelBuilder(ClientData clientData, Tcl_Interp *interp, int argc,
                        TCL_Char **argv);

extern int OPS_ResetInput(ClientData, Tcl_Interp *, int, int, TCL_Char **,
    Domain *, TclBuilder *);

int
myCommands(Tcl_Interp *interp)
{
  Tcl_CreateCommand(interp, "model", specifyModelBuilder, (ClientData)NULL,
                    (Tcl_CmdDeleteProc *)NULL);
  Tcl_Eval(interp, "rename load import;");
  return 0;
}

int
specifyModelBuilder(ClientData clientData, Tcl_Interp *interp, int argc,
                    TCL_Char **argv)
{
  /* int cArg = 0; */
  G3_Runtime *rt = G3_getRuntime(interp);
  TclBuilder *theNewBuilder = 0;
  Domain *theNewDomain = new Domain();
  G3_setDomain(rt,theNewDomain);

  // make sure at least one other argument to contain model builder type given
  if (argc < 2) {
    opserr << "WARNING need to specify a model type, valid types:\n";
    opserr << "\tBasicBuilder\n";
    return TCL_ERROR;
  }

  // invoke the descructor on the old builder
  if (theBuilder != 0) {
    delete theBuilder;
    theBuilder = 0;
  }

  // check argv[1] for type of ModelBuilder and create the object
  if ((strcmp(argv[1], "basic") == 0)        ||
      (strcmp(argv[1], "Basic") == 0)        ||
      (strcmp(argv[1], "safe") == 0)         ||
      (strcmp(argv[1], "BasicBuilder") == 0) ||
      (strcmp(argv[1], "basicBuilder") == 0)) {
    int ndm = 0;
    int ndf = 0;
    bool safe_builder = true;

    if (strcmp(argv[1], "safe") == 0) {
      safe_builder = true;
    }

    if (argc < 4) {
      opserr << "WARNING incorrect number of command arguments\n";
      opserr << "model modelBuilderType -ndm ndm? <-ndf ndf?> \n";
      return TCL_ERROR;
    }

    int posArg = 1; // track positional argument
    int argPos = 2;
    while (argPos < argc) {
      if (strcmp(argv[argPos], "-ndm") == 0 ||
          strcmp(argv[argPos], "-NDM") == 0) {
        argPos++;
        if (argPos < argc) {
          if (Tcl_GetInt(interp, argv[argPos], &ndm) != TCL_OK) {
            opserr << "WARNING error reading ndm: " << argv[argPos];
            opserr << "\nmodel modelBuilderType -ndm ndm? <-ndf ndf?>\n";
            return TCL_ERROR;
          }
        }
        argPos++;
        posArg++;
      } else if (strcmp(argv[argPos], "-ndf") == 0 ||
                 strcmp(argv[argPos], "-NDF") == 0) {
        argPos++;
        if (argPos < argc)
          if (Tcl_GetInt(interp, argv[argPos], &ndf) != TCL_OK) {
            opserr << "WARNING error reading ndf: " << argv[argPos];
            opserr << "\nmodel modelBuilderType -ndm ndm? <-ndf ndf?>\n";
            return TCL_ERROR;
          }
        argPos++;
        posArg++;
      } else if (posArg == 1) {
          if (Tcl_GetInt(interp, argv[argPos], &ndm) != TCL_OK) {
            opserr << "WARNING error reading ndm: " << argv[argPos];
            opserr << "\nmodel modelBuilderType -ndm ndm? <-ndf ndf?>\n";
            return TCL_ERROR;
          }
        argPos++;
        posArg++;
      } else if (posArg == 2) {
          if (Tcl_GetInt(interp, argv[argPos], &ndf) != TCL_OK) {
            opserr << "WARNING error reading ndf: " << argv[argPos];
            opserr << "\nmodel modelBuilderType -ndm ndm? <-ndf ndf?>\n";
            return TCL_ERROR;
          }
        argPos++;
        posArg++;
      } else {// Advance to next input argument if there are no matches -- MHS
        argPos++;
      }
    }

    // check that ndm was specified
    if (ndm == 0) {
      opserr << "WARNING need to specify ndm\n";
      opserr << "model modelBuilderType -ndm ndm? <-ndf ndf?>\n";
      return TCL_ERROR;
    }

    // check for ndf, if not assume one
    if (ndf == 0) {
      if (ndm == 1)
        ndf = 1;
      else if (ndm == 2)
        ndf = 3;
      else if (ndm == 3)
        ndf = 6;
      else {
        opserr << "WARNING specified ndm, " << ndm << ", will not work\n";
        opserr << "with any elements in BasicBuilder\n";
        return TCL_ERROR;
      }
    }
    // create the model builder
    theNewBuilder = new TclSafeBuilder(*theNewDomain, interp, ndm, ndf);

    if (theNewBuilder == 0) {
      opserr << "WARNING ran out of memory in creating BasicBuilder model\n";
      return TCL_ERROR;
    } else {
      theBuilder = theNewBuilder;
      G3_setModelBuilder(rt, theNewBuilder);
    }
  }
  else if ((strcmp(argv[1], "test") == 0) ||
           (strcmp(argv[1], "uniaxial") == 0) ||
           (strcmp(argv[1], "TestUniaxial") == 0) ||
           (strcmp(argv[1], "testUniaxial") == 0) ||
           (strcmp(argv[1], "UniaxialMaterialTest") == 0)) {
    int count = 1;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &count) != TCL_OK) {
        return TCL_ERROR;
      }
    }
    theNewBuilder = new TclUniaxialMaterialTester(theDomain, interp, count);
    if (theNewBuilder == 0) {
      opserr << "WARNING ran out of memory in creating "
                "TclUniaxialMaterialTester model\n";
      return TCL_ERROR;
    } else {
      G3_setModelBuilder(rt, theNewBuilder);
    }
  }
/*
  else if ((strcmp(argv[1], "testPlaneStress") == 0) ||
           (strcmp(argv[1], "PlaneStressMaterialTest") == 0)) {
    int count = 1;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &count) != TCL_OK) {
        return TCL_ERROR;
      }
    }

    theNewBuilder = new TclPlaneStressMaterialTester(theDomain, interp, count);
    if (theNewBuilder == 0) {
      opserr << "WARNING ran out of memory in creating "
                "TclUniaxialMAterialTester model\n";
      return TCL_ERROR;
    }
  }

*/
  else if ((strcmp(argv[1], "sectionTest") == 0) ||
           (strcmp(argv[1], "TestSection") == 0) ||
           (strcmp(argv[1], "testSection") == 0) ||
           (strcmp(argv[1], "SectionForceDeformationTest") == 0)) {
    int count = 1;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &count) != TCL_OK) {
        return TCL_ERROR;
      }
    }
    theNewBuilder = new TclSectionTestBuilder(theDomain, interp, count);
    if (theNewBuilder == 0) {
      opserr << "WARNING ran out of memory in creating "
                "TclUniaxialMAterialTester model\n";
      return TCL_ERROR;
    } else {
      G3_setModelBuilder(rt, theNewBuilder);
    }
  }
  else {
    opserr << "WARNING unknown model builder type\n";
    opserr << "WARNING model builder type " << argv[1] << " not supported\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}
