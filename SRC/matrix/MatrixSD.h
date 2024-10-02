//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//    See https://opensees.berkeley.edu/OpenSees/copyright.php for license.
//
//===--- MatrixSD.h - Matrix with fixed size ------------------------------===//

#pragma once

namespace OpenSees {

template <int n, bool half=false>
struct MatrixSD {
  constexpr static int nd = n*(n+1)/2;
  constexpr static int nr = n;

public:
  VectorND<n*(n+1)/2> vector;


  constexpr void
  zero()
  {
    for (int i = 0; i<nd; i++)
      vector[i] = 0.0;
  }

  constexpr MatrixND<n,n>
  full() const {
    MatrixND<n,n> S;
    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
        S(i,j) = (*this)(i,j);
    return S;
  }

  constexpr inline double& operator()(int i, int j) {
    return vector[vector_index(i, j, n)]; // *((half && i!=j)*(-0.5) + 1.0);
  }

  constexpr inline const double& operator()(int i, int j) const {
    return vector[vector_index(i, j, n)]*((half && i!=j)*(-0.5) + 1.0);
  }

  constexpr friend inline VectorND<n>
  operator*(const MatrixSD<n,half> &left, const VectorND<n> &right) {
    VectorND<n> prod;
    for (index_t i = 0; i < n; ++i) {
        prod[i] = 0.0;
        for (index_t k = 0; k < n; ++k) {
          prod[i] += left(i,k) * right[k];
        }
    }
    return prod;
  }


  template <index_t NR>
  inline constexpr friend MatrixND<NR, n>
  operator*(const MatrixND<NR, n> &left, const MatrixSD<n,half> &right) {
    MatrixND<NR, n> prod;
    for (index_t i = 0; i < NR; ++i) {
      for (index_t j = 0; j < n; ++j) {
        prod(i, j) = 0.0;
        for (index_t k = 0; k < n; ++k) {
          prod(i, j) += left(i,k) * right(k,j);
        }
      }
    }
    return prod;
  }

  template <index_t J>
  inline constexpr friend MatrixND<n, J>
  operator*(const MatrixSD<n,half> &left, const MatrixND<n, J> &right) {
    MatrixND<n, J> prod;
    for (index_t i = 0; i < n; ++i) {
      for (index_t j = 0; j < J; ++j) {
        prod(i, j) = 0.0;
        for (index_t k = 0; k < n; ++k) {
          prod(i, j) += left(i,k) * right(k,j);
        }
      }
    }
    return prod;
  }


  template <class MatT, int nk> void 
    addMatrixTransposeProduct(double thisFact, const MatrixND<nk, n, double> &, const MatT&, double scale);

  MatrixSD<n,half>& addDiagonal(const double vol);

private:

  constexpr static inline int
  vector_index(int i, int j, int N)
  {
     if (i <= j)
        return  i * N - (i - 1) * i / 2 + j - i;
     else
        return  j * N - (j - 1) * j / 2 + i - j;
  }
};


template <index_t n, bool h>
inline
MatrixSD<n,h>& 
MatrixSD<n,h>::addDiagonal(const double diag)
{
  for (int i=0; i<n; i++)
    vector[i] += diag;

  return *this;
}


// B'*C
template <int n, bool h>
template <class MatT, int nk> inline
void
MatrixSD<n,h>::addMatrixTransposeProduct(double thisFact,
                                       const MatrixND<nk, n>& B,
                                       const MatT& C,
                                       const double otherFact)
{
  MatrixSD<n,h>& self = *this;

  if (thisFact == 1.0) {
    for (int j=0; j<n; j++) {
      for (int i=0; i<n; i++) {
        const double *bkiPtr  = &(&B(0,0))[i*nk];
        const double *cjkPtr  = &(&C(0,0))[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        self(i,j) = self(i,j) + sum * otherFact;
      }
    }

  }

  // Set form
  else if (thisFact == 0.0) {
    for (int j=0; j<nr; j++) {
      for (int i=0; i<nr; i++) {
        const double *bkiPtr  = &(&B(0,0))[i*nk];
        const double *cjkPtr  = &(&C(0,0))[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        self(i,j) = sum * otherFact;
      }
    } 

  }

  // General form with BOTH thisFact and otherFact
  else {
    for (int j=0; j<nr; j++) {
      for (int i=0; i<nr; i++) {
        const double *bkiPtr  = &(&B(0,0))[i*nk];
        const double *cjkPtr  = &(&C(0,0))[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        self(i,j) = self(i,j) * thisFact + sum * otherFact;
      }
    } 
  }
}

} // namespace

