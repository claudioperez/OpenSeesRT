#include <array>
#include <string>
#include <unordered_map>

// framework
#include <g3_api.h>
// #include <runtimeAPI.h>
// #include <analysisAPI.h>
// #include <OPS_Globals.h>

// system of eqn and solvers
#include <BandSPDLinSOE.h>
#include <BandSPDLinLapackSolver.h>

#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>

#include <ConjugateGradientSolver.h>

#include <FullGenLinSOE.h>
#include <FullGenLinLapackSolver.h>

#include <ProfileSPDLinSOE.h>
#include <ProfileSPDLinDirectSolver.h>
#include <DistributedProfileSPDLinSOE.h>

#include <DiagonalSOE.h>
#include <DiagonalDirectSolver.h>

#include <SProfileSPDLinSolver.h>
#include <SProfileSPDLinSOE.h>

#include <SparseGenColLinSOE.h>

#include <SparseGenRowLinSOE.h>
// #include <SymSparseLinSOE.h>
// #include <SymSparseLinSolver.h>
#include <EigenSOE.h>
#include <EigenSolver.h>
#include <ArpackSOE.h>
#include <ArpackSolver.h>
#include <SymArpackSOE.h>
#include <SymArpackSolver.h>
#include <BandArpackSOE.h>
#include <BandArpackSolver.h>
#include <SymBandEigenSOE.h>
#include <SymBandEigenSolver.h>
#include <FullGenEigenSOE.h>
#include <FullGenEigenSolver.h>

#ifdef _THREADS
#  include <ThreadedSuperLU.h>
#else
#  include <SuperLU.h>
#endif

#ifdef _CUSP
#  include <CuSPSolver.h>
#endif

#ifdef _CULAS4
#include <CulaSparseSolverS4.h>
#endif

#ifdef _CULAS5
#include <CulaSparseSolverS5.h>
#endif

#ifdef _MUMPS
#  include <MumpsSOE.h>
#  include <MumpsSolver.h>
#  if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
#    include <MumpsParallelSOE.h>
#    include <MumpsParallelSolver.h>
#  endif
#endif

#if 1 || defined(_PETSC)
LinearSOE *TclCommand_newPetscSOE(int, TCL_Char**);
#endif

#ifdef _CUDA
#  include <BandGenLinSOE_Single.h>
#  include <BandGenLinLapackSolver_Single.h>
#endif

#if defined(_PARALLEL_PROCESSING)
//  parallel soe & solvers
#  include <DistributedBandSPDLinSOE.h>
#  include <DistributedSparseGenColLinSOE.h>
#  include <DistributedSparseGenRowLinSOE.h>
#  include <DistributedBandGenLinSOE.h>
#  include <DistributedDiagonalSOE.h>
#  include <DistributedDiagonalSolver.h>

#  include <MPIDiagonalSOE.h>
#  include <MPIDiagonalSolver.h>

#  define MPIPP_H
#  include <DistributedSuperLU.h>
#  include <DistributedProfileSPDLinSOE.h>
#endif

typedef LinearSOE* (fn)(G3_Runtime*, int, G3_Char**);

// Specifiers defined in solver.cpp
fn specify_SparseSPD;
fn specifySparseGen;
fn G3Parse_newMumpsLinearSOE;


// Helpers to automatically create constructors for systems/solvers 
// that do not take arguments when they are constructed.
template <typename Solver, typename SOE>
LinearSOE *simple_soe(G3_Runtime*, int, G3_Char**) {return new SOE(*(new Solver()));}

#define MK_SOE(Solver, SOE) simple_soe<Solver, SOE>

std::unordered_map<std::string, fn*> LinearSOE_Library = {
  {"BandSPD",               MK_SOE(BandSPDLinLapackSolver,      BandSPDLinSOE)},

  {"BandGeneral",           MK_SOE(BandGenLinLapackSolver,      BandGenLinSOE)},
  // BandGen, BandGEN

  {"BandSPD",               MK_SOE(BandSPDLinLapackSolver,      BandSPDLinSOE)},

  {"SparseGen",             specifySparseGen},
  {"SuperLU",               specifySparseGen},

  {"SparseSPD",             specify_SparseSPD},

  {"Diagonal",              MK_SOE(DiagonalDirectSolver,        DiagonalSOE)},

  {"MPIDiagonal",           MK_SOE(DiagonalDirectSolver,        DiagonalSOE)},

  {"SProfileSPD",           MK_SOE(SProfileSPDLinSolver,        SProfileSPDLinSOE)},

  {"ProfileSPD",            MK_SOE(ProfileSPDLinDirectSolver,   ProfileSPDLinSOE)},

  {"FullGeneral",           MK_SOE(FullGenLinLapackSolver,      FullGenLinSOE)},

#ifdef _CUDA
  {"BandGeneral_Single", MK_SOE(BandGenLinLapackSolver_Single,    BandGenLinSOE_Single)},
#endif
};


