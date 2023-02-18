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
//
// Description: This file contains the class interface for ThreadedSuperLU.
// This is a class that uses the threads version of SuperLU
//
// File: ~/system_of_eqn/linearSOE/sparseGEN/ThreadedSuperLU.h
//
// Written: fmk 
// Created: 11/96
//
#ifndef ThreadedSuperLU_h
#define ThreadedSuperLU_h

#define v3

#ifdef v3
#  include "SuperLU_MT-3/SRC/slu_mt_ddefs.h"
#  include "SuperLU_MT-3/SRC/supermatrix.h"
#  define pdgstrf_options_t superlumt_options_t
#  define NC SLU_NC
#  define GE SLU_GE
#  define DN SLU_DN
#  define _D SLU_D
#else
// #include <pdsp_defs.h>
#  include "SuperLU_MT/pdsp_defs.h"
#  include "SuperLU_MT/supermatrix.h"
#endif

#include <SparseGenColLinSolver.h>

class ThreadedSuperLU : public SparseGenColLinSolver
{
  public:
    ThreadedSuperLU(int numThreads = 2,
		    int permSpec = 0, int panelSize = 6, 
		    int relax = 6, double thresh = 0.0);     

    ~ThreadedSuperLU();

    int solve(void);
    int setSize(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);        
    
  protected:

  private:
    SuperMatrix A,L,U,B,AC;
    int *perm_r;
    int *perm_c;
    int *etree;
    int sizePerm;
    int relax, permSpec, panelSize;
    float thresh;

    int numThreads;
    pdgstrf_options_t pdgstrf_options;
    yes_no_t refact, usepr;
    fact_t refact_3; // for version 3
    fact_t fact;
    trans_t trans;
    void *work;
    int lwork;
    Gstat_t gStat;
};

#endif


