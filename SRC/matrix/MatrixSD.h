#pragma once

namespace OpenSees {

template <int n>
struct MatrixSD {
  constexpr static int nd = n*(n+1)/2;
  constexpr static int nr = n;

  VectorND<n*(n+1)/2> vector;

public:

  constexpr double& operator()(int i, int j) {
    return vector[vector_index(i, j, n*(n+1)/2)];
  }

  constexpr const double& operator()(int i, int j) const {
    return vector[vector_index(i, j, n*(n+1)/2)];
  }



  template <class MatT, int nk> void 
    addMatrixTransposeProduct(double thisFact, const MatrixND<nk, n, double> &, const MatT&, double scale);

  MatrixSD<n>& addDiagonal(const double vol);

private:

  constexpr static inline int
  vector_index(int i, int j, int N)
  {
     if (i <= j)
        return i * N - (i - 1) * i / 2 + j - i;
     else
        return j * N - (j - 1) * j / 2 + i - j;
  }
};


template <index_t n> inline
MatrixSD<n>& 
MatrixSD<n>::addDiagonal(const double diag)
{
  for (int i=0; i<n; i++)
    (*this)(i,i) += diag;

  return *this;
}

// B'*C
template <int n> 
template <class MatT, int nk> inline
void
MatrixSD<n>::addMatrixTransposeProduct(double thisFact,
                                       const MatrixND<nk, n>& B,
                                       const MatT& C,
                                       const double otherFact)
{
  if (thisFact == 1.0) {
    double *aijPtr = &(*this)(0,0);
    for (int j=0; j<n; j++) {
      for (int i=0; i<n; i++) {
        const double *bkiPtr  = &(&B(0,0))[i*nk];
        const double *cjkPtr  = &(&C(0,0))[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        *aijPtr++ += sum * otherFact;
      }
    }

  } else if (thisFact == 0.0) {
    double *aijPtr = &(*this)(0,0);
    for (int j=0; j<nr; j++) {
      for (int i=0; i<nr; i++) {
        const double *bkiPtr  = &(&B(0,0))[i*nk];
        const double *cjkPtr  = &(&C(0,0))[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        *aijPtr++ = sum * otherFact;
      }
    } 
  } else {
    double *aijPtr = &(*this)(0,0);
    for (int j=0; j<nr; j++) {
      for (int i=0; i<nr; i++) {
        const double *bkiPtr  = &(&B(0,0))[i*nk];
        const double *cjkPtr  = &(&C(0,0))[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        *aijPtr = *aijPtr * thisFact + sum * otherFact;
        aijPtr++;
      }
    } 
  }
}

} // namespace

