//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//    See https://opensees.berkeley.edu/OpenSees/copyright.php for license.
//
//===--- MatrixSD.h - Matrix with fixed size ------------------------------===//

#pragma once
#include <cassert>

namespace OpenSees {

template <int n, bool half=false>
struct MatrixSD {
  constexpr static int nd = n*(n+1)/2;
  constexpr static int nr = n;
  constexpr static int size = n*(n+1)/2;

  VectorND<n*(n+1)/2> vector;


  constexpr inline void zero();
  constexpr inline MatrixND<n,n> full() const;

//constexpr inline double& operator()(int i, int j) {
//  return vector[vector_index(i, j)]; // *((half && i!=j)*(-0.5) + 1.0);
//}

  constexpr inline double& ref(int i, int j) {
    return vector[vector_index(i, j)]; // *((half && i!=j)*(-0.5) + 1.0);
  }

  // Getter for Voigt notation with division for off-diagonal elements
  constexpr inline double  operator()(int i, int j) const {
    return vector[vector_index(i, j)]*((half && i!=j)*(-0.5) + 1.0);
  }

  template <class MatT, int nk> void 
    addMatrixTransposeProduct(double thisFact, 
                              const MatrixND<nk, n, double> &, 
                              const MatT&, double scale);

  constexpr inline MatrixSD<n,half>& addDiagonal(double vol);

  template <typename MatTyp>
  constexpr MatrixSD<n,half> &
  operator=(const MatTyp &other)
  {
    for (index_t j = 0; j < n; ++j) {
      for (index_t i = 0; i < n; ++i) {
        ref(i,j) = other(i,j);
      }
    }
    return *this;
  }

  constexpr MatrixND<n,half> &
  operator+=(const double value) {
    vector += value;
    return *this;
  }

  constexpr MatrixND<n,half> &
  operator*=(const double value) {
    vector *= value;
    return *this;
  }

  // FRIENDS
#if 0
  friend constexpr MatrixSD
  operator+(MatrixSD left, const MatrixSD &right) {
    left += right; 
    return left;
  }

  friend constexpr MatrixSD
  operator-(MatrixSD left, const MatrixSD &right) {
    left -= right; 
    return left;
  }
#endif
  
  friend constexpr MatrixSD<n,half> // scalar * Matrix
  operator*(double scalar, MatrixSD<n,half> mat) {
    mat *= scalar;
    return mat;
  }

  friend constexpr MatrixSD<n,half> // Matrix * scalar
  operator*(MatrixSD<n,half> mat, double scalar) {
    mat *= scalar;
    return mat;
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


  // Function to get the index in the vector for A(i, j)
  constexpr static inline int
  vector_index(int i, int j)
  {
    assert(i >= 0 && i < n && j >= 0 && j < n);

    if (i == j) {
        return i; // Diagonal elements
    } else if (i < j) {
        // Upper triangle
        return n + (i * (2 * n - i - 1) / 2 + (j - i - 1)); 
    } else {
        // Lower triangle
        return n + (j * (2 * n - j - 1) / 2 + (i - j - 1)); 
    }
//   if (i <= j)
//      return  i * N - (i - 1) * i / 2 + j - i;
//   else
//      return  j * N - (j - 1) * j / 2 + i - j;
  }
};

//
//
//

template <int n, bool h>
constexpr void
MatrixSD<n,h>::zero()
{
  for (int i = 0; i<nd; i++)
    vector[i] = 0.0;
}


template <int n, bool h>
constexpr MatrixND<n,n>
MatrixSD<n,h>::full() const 
{
    MatrixND<n,n> S;
    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
        S(i,j) = (*this)(i,j);
    return S;
  }


template <index_t n, bool h>
constexpr MatrixSD<n,h>& 
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
//static_assert(!h, "Method not supported for Voigt type");

  MatrixSD<n,h>& self = *this;

  if (thisFact == 1.0) {
    for (int j=0; j<n; j++) {
      for (int i=0; i<=j; i++) {
        const double *bkiPtr  = &(&B(0,0))[i*nk];
        const double *cjkPtr  = &(&C(0,0))[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        ref(i,j) = self(i,j) + sum * otherFact;
      }
    }

  }

  // Set form
  else if (thisFact == 0.0) {
    for (int j=0; j<n; j++) {
      for (int i=0; i<=j; i++) {
        const double *bkiPtr  = &(&B(0,0))[i*nk];
        const double *cjkPtr  = &(&C(0,0))[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        ref(i,j) = sum * otherFact;
      }
    } 

  }

  // General form with BOTH thisFact and otherFact
  else {
    for (int j=0; j<n; j++) {
      for (int i=0; i<=j; i++) {
        const double *bkiPtr  = &(&B(0,0))[i*nk];
        const double *cjkPtr  = &(&C(0,0))[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        ref(i,j) = self(i,j) * thisFact + sum * otherFact;
      }
    } 
  }
}

} // namespace

