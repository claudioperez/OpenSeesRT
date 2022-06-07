#include <G3Parse.h>
#include <ItpackLinSOE.h>
#include <ItpackLinSolver.h>

LinearSOE*
G3Parse_newItpackLinearSOE(G3_Builder *rt, int argc, TCL_Char **argv)
{
    // now must determine the type of solver to create 
    // from rest of args
    int method = 1;
    if (argc == 3) {
      if (G3Parse_getInt(rt, argv[2], &method) != TCL_OK)
        return nullptr;
    }
    ItpackLinSolver *theSolver = new ItpackLinSolver(method);
    return  new ItpackLinSOE(*theSolver);
}



