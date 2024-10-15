//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the implementation for the FeapElmt class.
//
// Written: cmp
// 
#include "fElement.h"
#include <ID.h>
#include <Vector.h>

typedef int elmt_t(double *d, double *ul, double *xl, int *ix, double *tl, 
                   double *s, double *r, int *ndf, int *ndm, int *nst, int *isw, 
                   double *dm, int *nen, int *n, int *nh1, int *nh2, int *nh3, 
                   double *h, double *ctan, int *ior, int *iow);


template <elmt_t *elmt, const char* name, int NST, int NEN, int NDF, int NDM>
class FeapElmt : public fElement {
private:

public:

  FeapElmt(int tag, int nn, int *nodes, int nd, double* d)
    : fElement(tag, ELE_TAG_FeapElmt, 5, nd, NEN, NDM, NDF, 0, 0)
  {
    for (int i=0; i<nd; i++)
    (*data)(i) = d[i];

    for (int i=0; i<nn; i++)
      (*connectedNodes)(i) = nodes[i];
  }
      
  FeapElmt()
  :fElement(ELE_TAG_FeapElmt)    
  {
      // does nothing
  }

  ~FeapElmt()
  {
      // does nothing
  }

  bool 
  validate()
  {
    return true;
  }

  int
  invokefRoutine(double *d, double *ul, double *xl, int *ix, double *tl, 
                 double *s, double *r, int ndf, int ndm, int nst, int isw, 
                 double dm, int nen, int n, int nh1, int nh2, int nh3, 
                 double *h, double *ctan, int ior, int iow)
  {
      // check that the values are acceptable to the fortran subroutine
      if (nst != NST || nen != NEN || dm != 2)
          return 0;
      
      elmt(d, ul, xl, ix, tl, s, r, &ndf, &ndm, &nst, &isw, &dm,
              &nen, &n, &nh1, &nh2, &nh3, h, ctan, &ior, &iow);
          
      return nst;
  }
};

