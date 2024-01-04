/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Claudio M. Perez
//
#include <stdio.h>
#include <stdlib.h>
#include "GaussLegendre1D.hpp"

template<template<int, int> class Gauss> void
gauss(int n, const double **const x, const double **const w)
{
  switch (n) {
    case  1: *x = Gauss<1, 1>::pts; *w = Gauss<1, 1>::wts; break;
    case  2: *x = Gauss<1, 2>::pts; *w = Gauss<1, 2>::wts; break;
    case  3: *x = Gauss<1, 3>::pts; *w = Gauss<1, 3>::wts; break;
    case  4: *x = Gauss<1, 4>::pts; *w = Gauss<1, 4>::wts; break;
    case  5: *x = Gauss<1, 5>::pts; *w = Gauss<1, 5>::wts; break;
    case  6: *x = Gauss<1, 6>::pts; *w = Gauss<1, 6>::wts; break;
    case  7: *x = Gauss<1, 7>::pts; *w = Gauss<1, 7>::wts; break;
    case  8: *x = Gauss<1, 8>::pts; *w = Gauss<1, 8>::wts; break;
    case  9: *x = Gauss<1, 9>::pts; *w = Gauss<1, 9>::wts; break;
    case 10: *x = Gauss<1,10>::pts; *w = Gauss<1,10>::wts; break;
    case 11: *x = Gauss<1,11>::pts; *w = Gauss<1,11>::wts; break;
    case 12: *x = Gauss<1,12>::pts; *w = Gauss<1,12>::wts; break;
    case 13: *x = Gauss<1,13>::pts; *w = Gauss<1,13>::wts; break;
    case 14: *x = Gauss<1,14>::pts; *w = Gauss<1,14>::wts; break;
    case 15: *x = Gauss<1,15>::pts; *w = Gauss<1,15>::wts; break;
    case 16: *x = Gauss<1,16>::pts; *w = Gauss<1,16>::wts; break;
    case 17: *x = Gauss<1,17>::pts; *w = Gauss<1,17>::wts; break;
    case 18: *x = Gauss<1,18>::pts; *w = Gauss<1,18>::wts; break;
    case 19: *x = Gauss<1,19>::pts; *w = Gauss<1,19>::wts; break;
    case 20: *x = Gauss<1,20>::pts; *w = Gauss<1,20>::wts; break;
  }
}

int
main(int argc, char** argv)
{
  const double *x, *w;
  int n = atoi(argv[1]);

  gauss<GaussLegendre>(n, &x, &w);

  for (int i=0; i<n; i++)
    printf("%lf\t%lf\n", x[i], w[i]);
  
}

