#include <tcl.h>
#include <string.h>
#include <UmfpackGenLinSOE.h>
#include <UmfpackGenLinSolver.h>


LinearSOE*
TclDispatch_newUmfpackLinearSOE(ClientData clientData, Tcl_Interp* interp, int argc, const char** const argv)
{
    // now must determine the type of solver to create 
    // from rest of args
    int factLVALUE = 10;
    int factorOnce = 0;
    int printTime = 0;
    int count = 2;

    while (count < argc) {
      if ((strcmp(argv[count], "-lValueFact") == 0) ||
          (strcmp(argv[count], "-lvalueFact") == 0) ||
          (strcmp(argv[count], "-LVALUE") == 0)) {
        if (count+1 < argc && Tcl_GetInt(interp, argv[count + 1], &factLVALUE) != TCL_OK)
          return nullptr;
        count++;
      } else if ((strcmp(argv[count], "-factorOnce") == 0) ||
                 (strcmp(argv[count], "-FactorOnce") == 0)) {
        factorOnce = 1;
        count++;
      } else if ((strcmp(argv[count], "-printTime") == 0) ||
                 (strcmp(argv[count], "-time") == 0)) {
        printTime = 1;
        count++;
      }
    }
    UmfpackGenLinSolver *theSolver = new UmfpackGenLinSolver();
    // theSOE = new UmfpackGenLinSOE(*theSolver, factLVALUE, factorOnce, printTime);
    return new UmfpackGenLinSOE(*theSolver);
}

