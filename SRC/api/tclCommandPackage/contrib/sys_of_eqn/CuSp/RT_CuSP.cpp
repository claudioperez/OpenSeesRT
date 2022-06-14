
#ifdef _CUSP
#  include <CuSPSolver.h>
#endif


#ifdef _CUSP
else if ((_stricmp(argv[1], "CuSP") == 0)) {

  double relTol = 1e-6;
  int maxInteration = 100000;
  int preCond = 1; // diagonal
  int solver = 0;  // Bicg
  int count = 2;

  while (count < argc) {
    if (_stricmp(argv[count], "-rTol") == 0) {
      count++;
      if (count < argc)
        if (Tcl_GetDouble(interp, argv[count], &relTol) != TCL_OK)
          return TCL_ERROR;
    } else if ((_stricmp(argv[count], "-mInt") == 0)) {
      count++;
      if (count < argc)
        if (Tcl_GetInt(interp, argv[count], &maxInteration) != TCL_OK)
          return TCL_ERROR;
    } else if ((_stricmp(argv[count], "-pre") == 0)) {
      count++;
      if (count < argc)
        if ((_stricmp(argv[count], "none") == 0))
          preCond = 0;
        else if ((_stricmp(argv[count], "diagonal") == 0))
          preCond = 1;
        else if ((_stricmp(argv[count], "ainv") == 0))
          preCond = 2;
        else
          return TCL_ERROR;
    } else if ((_stricmp(argv[count], "-solver") == 0)) {
      count++;
      if (count < argc)
        if ((_stricmp(argv[count], "bicg") == 0))
          solver = 0;
        else if ((_stricmp(argv[count], "bicgstab") == 0))
          solver = 1;
        else if ((_stricmp(argv[count], "cg") == 0))
          solver = 2;
        else if ((_stricmp(argv[count], "gmres") == 0))
          solver = 3;
        else
          return TCL_ERROR;
    }
    count++;
  }

  CuSPSolver *theSolver =
      new CuSPSolver(maxInteration, relTol, preCond, solver);
  theSOE = new SparseGenRowLinSOE(*theSolver);
}
#endif // _CUSP
