#include <ItpackLinSOE.h>
#include <ItpackLinSolver.h>

TclCommand_addItpack(int argc, TCL_Char **argv)
{
    // now must determine the type of solver to create 
    // from rest of args
    int method = 1;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &method) != TCL_OK)
        return TCL_ERROR;
    }
    ItpackLinSolver *theSolver = new ItpackLinSolver(method);
    theSOE = new ItpackLinSOE(*theSolver);
  }
}



