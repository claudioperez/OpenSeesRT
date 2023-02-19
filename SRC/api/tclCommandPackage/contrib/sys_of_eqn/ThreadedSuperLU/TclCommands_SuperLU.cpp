
#include <tcl.h>
#ifdef _THREADS
#  include "contrib/sys_of_eqn/SuperLU/ThreadedSuperLU.h"
#else
#endif
// TODO: CMP
// #include <SuperLU.h>

LinearSOE*
specifySparseGen(G3_Runtime* rt, int argc, G3_Char **argv)
{
  // SPARSE GENERAL SOE * SOLVER
  if ((strcmp(argv[1], "SparseGeneral") == 0) ||
           (strcmp(argv[1], "SuperLU") == 0) ||
           (strcmp(argv[1], "SparseGEN") == 0)) {
    Tcl_Interp *interp = G3_getInterpreter(rt);

    SparseGenColLinSolver *theSolver = 0;
    int count = 2;
    double thresh = 0.0;
    int npRow = 1;
    int npCol = 1;
    int np = 1;

    // defaults for threaded SuperLU
    while (count < argc) {
      if ((strcmp(argv[count], "p")    == 0) ||
          (strcmp(argv[count], "piv")  == 0) ||
          (strcmp(argv[count], "-piv") == 0)) {
        thresh = 1.0;
      } else if ((strcmp(argv[count], "-np") == 0) ||
                 (strcmp(argv[count], "np")  == 0)) {
        count++;
        if (count < argc)
          if (Tcl_GetInt(interp, argv[count], &np) != TCL_OK)
            return nullptr;
      } else if ((strcmp(argv[count], "npRow")  == 0) ||
                 (strcmp(argv[count], "-npRow") == 0)) {
        count++;
        if (count < argc)
          if (Tcl_GetInt(interp, argv[count], &npRow) != TCL_OK)
            return nullptr;
      } else if ((strcmp(argv[count], "npCol")  == 0) ||
                 (strcmp(argv[count], "-npCol") == 0)) {
        count++;
        if (count < argc)
          if (Tcl_GetInt(interp, argv[count], &npCol) != TCL_OK)
            return nullptr;
      }
      count++;
    }

    int permSpec = 0;
    int panelSize = 6;
    int relax = 6;

#ifdef _THREADS
    if (np != 0)
      theSolver = new ThreadedSuperLU(np, permSpec, panelSize, relax, thresh);
#endif

#ifdef _PARALLEL_PROCESSING
    if (theSolver != 0)
      delete theSolver;
    theSolver = 0;

    if (npRow != 0 && npCol != 0) {
      theSolver = new DistributedSuperLU(npRow, npCol);
    }
#else
    char symmetric = 'N';
    double drop_tol = 0.0;
    while (count < argc) {
      if (strcmp(argv[count], "s") == 0 || strcmp(argv[count], "symmetric") ||
          strcmp(argv[count], "-symm")) {
        symmetric = 'Y';
      }
      count++;
    }
    // TODO : CMP
    // theSolver = new SuperLU(permSpec, drop_tol, panelSize, relax, symmetric);
#endif

#ifdef _PARALLEL_PROCESSING
    return new DistributedSparseGenColLinSOE(*theSolver);
#else
    return new SparseGenColLinSOE(*theSolver);
#endif

  } else {
    return nullptr;
  }
}

