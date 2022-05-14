/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.1 $
// $Date: 2005-08-04 00:18:03 $
// $Source:
// /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/petsc/TclPetsc.cpp,v $

// Written: fmk

// Description: This file contains the function invoked when the user invokes
// the Petsc command in the interpreter.

// #include <OPS_Globals.h>
#include <stdlib.h>
#include <string.h>
#include <g3_api.h>

#include <PetscSOE.h>
#include <PetscSolver.h>
#include <PetscSparseSeqSolver.h>
#include <SparseGenRowLinSOE.h>


#ifdef _USRDLL
#define DllExport _declspec(dllexport)
#else
#define DllExport
#endif

LinearSOE*
TclCommand_newPetscSOE(int argc, TCL_Char **argv) //, LinearSOE **theSOE)
{
  LinearSOE** theSOE;

  // now must determine the type of solver to create from rest of args
  KSPType method = KSPCG;           // KSPCG KSPGMRES
  PCType preconditioner = PCJACOBI; // PCJACOBI PCILU PCBJACOBI
  int matType = 0;

  double rTol = 1.0e-5;
  double aTol = 1.0e-50;
  double dTol = 1.0e5;
  int maxIts = 100000;
  int count = 2;
  while (count < argc - 1) {
    if (strcmp(argv[count], "-matrixType") == 0 ||
        strcmp(argv[count], "-matrix")) {
      if (strcmp(argv[count + 1], "sparse") == 0)
        matType = 1;
    } else if (strcmp(argv[count], "-rTol") == 0 ||
               strcmp(argv[count], "-relTol") ||
               strcmp(argv[count], "-relativeTolerance")) {
      if (Tcl_GetDouble(interp, argv[count + 1], &rTol) != TCL_OK)
        return nullptr;
    } else if (strcmp(argv[count], "-aTol") == 0 ||
               strcmp(argv[count], "-absTol") ||
               strcmp(argv[count], "-absoluteTolerance")) {
      if (Tcl_GetDouble(interp, argv[count + 1], &aTol) != TCL_OK)
        return nullptr;
    } else if (strcmp(argv[count], "-dTol") == 0 ||
               strcmp(argv[count], "-divTol") ||
               strcmp(argv[count], "-divergenceTolerance")) {
      if (Tcl_GetDouble(interp, argv[count + 1], &dTol) != TCL_OK)
        return nullptr;
    } else if (strcmp(argv[count], "-mIts") == 0 ||
               strcmp(argv[count], "-maxIts") ||
               strcmp(argv[count], "-maxIterations")) {
      if (Tcl_GetInt(interp, argv[count + 1], &maxIts) != TCL_OK)
        return nullptr;
    } else if (strcmp(argv[count], "-KSP") == 0 ||
               strcmp(argv[count], "-KSPType")) {
      if (strcmp(argv[count + 1], "KSPCG") == 0)
        method = KSPCG;
      else if (strcmp(argv[count + 1], "KSPBICG") == 0)
        method = KSPBICG;
      else if (strcmp(argv[count + 1], "KSPRICHARDSON") == 0)
        method = KSPRICHARDSON;
      else if (strcmp(argv[count + 1], "KSPCHEBYCHEV") == 0)
        method = KSPCHEBYCHEV;
      else if (strcmp(argv[count + 1], "KSPGMRES") == 0)
        method = KSPGMRES;
    } else if (strcmp(argv[count], "-PC") == 0 ||
               strcmp(argv[count], "-PCType")) {
      if ((strcmp(argv[count + 1], "PCJACOBI") == 0) ||
          (strcmp(argv[count + 1], "JACOBI") == 0))
        preconditioner = PCJACOBI;
      else if ((strcmp(argv[count + 1], "PCILU") == 0) ||
               (strcmp(argv[count + 1], "ILU") == 0))
        preconditioner = PCILU;
      else if ((strcmp(argv[count + 1], "PCICC") == 0) ||
               (strcmp(argv[count + 1], "ICC") == 0))
        preconditioner = PCICC;
      else if ((strcmp(argv[count + 1], "PCBJACOBI") == 0) ||
               (strcmp(argv[count + 1], "BIJACOBI") == 0))
        preconditioner = PCBJACOBI;
      else if ((strcmp(argv[count + 1], "PCNONE") == 0) ||
               (strcmp(argv[count + 1], "NONE") == 0))
        preconditioner = PCNONE;
    }
    count += 2;
  }

  if (matType == 0) {
    PetscSolver *theSolver =
        new PetscSolver(method, preconditioner, rTol, aTol, dTol, maxIts);
    (*theSOE) = new PetscSOE(*theSolver);
  } else {
    PetscSparseSeqSolver *theSolver = new PetscSparseSeqSolver(
        method, preconditioner, rTol, aTol, dTol, maxIts);
    (*theSOE) = new SparseGenRowLinSOE(*theSolver);
  }

  if (*theSOE == 0) {
    opserr << "WARNING system Petsc - out of memory\n";
    return nullptr;
  }

  return *theSOE;
}

#if defined(OPSDEF_USE_PETSC)
  else if (strcmp(argv[1], "Petsc") == 0) {
    // now must determine the type of solver to create from rest of args
    KSPType method = KSPCG;           // KSPCG KSPGMRES
    PCType preconditioner = PCJACOBI; // PCJACOBI PCILU PCBJACOBI
    int matType = 0;

    double rTol = 1.0e-5;
    double aTol = 1.0e-50;
    double dTol = 1.0e5;
    int maxIts = 100000;
    int count = 2;
    while (count < argc - 1) {
      if (strcmp(argv[count], "-matrixType") == 0 ||
          strcmp(argv[count], "-matrix")) {
        if (strcmp(argv[count + 1], "sparse") == 0)
          matType = 1;
      } else if (strcmp(argv[count], "-rTol") == 0 ||
                 strcmp(argv[count], "-relTol")    ||
                 strcmp(argv[count], "-relativeTolerance")) {
        if (Tcl_GetDouble(interp, argv[count + 1], &rTol) != TCL_OK)
          return nullptr;
      } else if (strcmp(argv[count], "-aTol") == 0 ||
                 strcmp(argv[count], "-absTol") ||
                 strcmp(argv[count], "-absoluteTolerance")) {
        if (Tcl_GetDouble(interp, argv[count + 1], &aTol) != TCL_OK)
          return nullptr;
      } else if (strcmp(argv[count], "-dTol") == 0 ||
                 strcmp(argv[count], "-divTol") ||
                 strcmp(argv[count], "-divergenceTolerance")) {
        if (Tcl_GetDouble(interp, argv[count + 1], &dTol) != TCL_OK)
          return nullptr;
      } else if (strcmp(argv[count], "-mIts") == 0 ||
                 strcmp(argv[count], "-maxIts") ||
                 strcmp(argv[count], "-maxIterations")) {
        if (Tcl_GetInt(interp, argv[count + 1], &maxIts) != TCL_OK)
          return nullptr;
      } else if (strcmp(argv[count], "-KSP") == 0 ||
                 strcmp(argv[count], "-KSPType")) {
        if (strcmp(argv[count + 1], "KSPCG") == 0)
          method = KSPCG;
        else if (strcmp(argv[count + 1], "KSPBICG") == 0)
          method = KSPBICG;
        else if (strcmp(argv[count + 1], "KSPRICHARDSON") == 0)
          method = KSPRICHARDSON;
        else if (strcmp(argv[count + 1], "KSPCHEBYSHEV") == 0)
          method = KSPCHEBYSHEV;
        else if (strcmp(argv[count + 1], "KSPGMRES") == 0)
          method = KSPGMRES;
      } else if (strcmp(argv[count], "-PC") == 0 ||
                 strcmp(argv[count], "-PCType")) {
        if ((strcmp(argv[count + 1], "PCJACOBI") == 0) ||
            (strcmp(argv[count + 1], "JACOBI") == 0))
          preconditioner = PCJACOBI;
        else if ((strcmp(argv[count + 1], "PCILU") == 0) ||
                 (strcmp(argv[count + 1], "ILU") == 0))
          preconditioner = PCILU;
        else if ((strcmp(argv[count + 1], "PCICC") == 0) ||
                 (strcmp(argv[count + 1], "ICC") == 0))
          preconditioner = PCICC;
        else if ((strcmp(argv[count + 1], "PCBJACOBI") == 0) ||
                 (strcmp(argv[count + 1], "BIJACOBI") == 0))
          preconditioner = PCBJACOBI;
        else if ((strcmp(argv[count + 1], "PCNONE") == 0) ||
                 (strcmp(argv[count + 1], "NONE") == 0))
          preconditioner = PCNONE;
      }
      count += 2;
    }

    if (matType == 0) {
      // PetscSolver *theSolver = new PetscSolver(method, preconditioner, rTol,
      // aTol, dTol, maxIts);
      PetscSolver *theSolver =
          new PetscSolver(method, preconditioner, rTol, aTol, dTol, maxIts);
      theSOE = new PetscSOE(*theSolver);
    } else {
      // PetscSparseSeqSolver *theSolver = new PetscSparseSeqSolver(method,
      // preconditioner, rTol, aTol, dTol, maxIts);
      PetscSparseSeqSolver *theSolver = 0;
      theSOE = new SparseGenRowLinSOE(*theSolver);
    }
  }
#endif
