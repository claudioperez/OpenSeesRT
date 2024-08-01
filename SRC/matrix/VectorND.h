//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// 
//  Objectives:
//  - little to no overhead above C-style arrays
//  - value semantics; objects do not decay to pointers;
//
//  This code is influenced by the following sources
//   list initialization:
//   - https://stackoverflow.com/questions/42068882/list-initialization-for-a-matrix-class
//
//   style/practices
//   - https://quuxplusone.github.io/blog/2021/04/03/static-constexpr-whittling-knife/
//  
//   Operator overloading / semantics
//   - https://stackoverflow.com/questions/9851188/does-it-make-sense-to-use-move-semantics-for-operator-and-or-operator/9851423#9851423
//
//   compile-time template restrictions/concepts:
//   - https://codereview.stackexchange.com/questions/259038/compile-time-matrix-class
//     (C++ 20)
//   - https://github.com/calebzulawski/cotila/
//     (C++ 17)
//
//===----------------------------------------------------------------------===//
//
// Claudio Perez
//
#ifndef VectorND_H
#define VectorND_H
#include <math.h>
#include <assert.h>
#include <array>
#include <stdexcept>
#include <functional>
#include <iostream> // overloading <<
#include <Vector.h>
#include <Matrix.h>
#include "blasdecl.h"

#if __cplusplus < 202000L
#  define consteval
#  define requires(X)
#endif

namespace OpenSees {

typedef int index_t;

template<int n, int m, typename T> struct MatrixND;

template <index_t N, typename T=double> 
requires(N > 0)
struct VectorND {
  T values[N];

  operator Vector() { return Vector(values, N);}

  template<int n, int m, typename> friend struct MatrixND;

  inline constexpr T&
  operator[](index_t index) {
    return values[index];
  }

  inline constexpr const T&
  operator[](index_t index) const {
    return values[index];
  }

  inline constexpr T&
  operator()(index_t index) {
    return values[index];
  }

  inline constexpr const T&
  operator()(index_t index) const {
    return values[index];
  }

  consteval int
  size() const {
    return N;
  }

  consteval inline void
  zero() {
    for (T& item : values )
      item = 0.0;
//  for (index_t i = 0; i < N; ++i)
//    values[i] = 0.0;
  }

  int
  addVector(const T thisFact, const Vector &other, const T otherFact) {
    if (otherFact == 0.0 && thisFact == 1.0)
      return 0; 

    else if (thisFact == 1.0) {
      // want: this += other * otherFact
      double *dataPtr = values;
      double *otherDataPtr = other.theData;
      if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) 
          *dataPtr++ += *otherDataPtr++;
      } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) 
          *dataPtr++ -= *otherDataPtr++;
      } else 
        for (int i=0; i<N; i++) 
          *dataPtr++ += *otherDataPtr++ * otherFact;

    } else if (thisFact == 0.0) {
        // want: this = other * otherFact
        double *dataPtr = values;
        double *otherDataPtr = other.theData;
        if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
          for (int i=0; i<N; i++) 
            *dataPtr++ = *otherDataPtr++;
        } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
          for (int i=0; i<N; i++) 
            *dataPtr++ = -(*otherDataPtr++);
        } else 
          for (int i=0; i<N; i++) 
            *dataPtr++ = *otherDataPtr++ * otherFact;
    } else {
      // want: this = this * thisFact + other * otherFact
      double *dataPtr = values;
      double *otherDataPtr = other.theData;
      if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) {
          double value = *dataPtr * thisFact + *otherDataPtr++;
          *dataPtr++ = value;
        }
      } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) {
          double value = *dataPtr * thisFact - *otherDataPtr++;
          *dataPtr++ = value;
        }
      } else 
        for (int i=0; i<N; i++) {
          double value = *dataPtr * thisFact + *otherDataPtr++ * otherFact;
          *dataPtr++ = value;
      }
    }

    // successfull
    return 0;
  }

  int
  addVector(const T thisFact, const VectorND<N> &other, const T otherFact) {
    if (otherFact == 0.0 && thisFact == 1.0)
      return 0; 

    else if (thisFact == 1.0) {
      // want: this += other * otherFact
      double *dataPtr = values;
      const double * otherDataPtr = other.values;
      if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) 
          *dataPtr++ += *otherDataPtr++;
      } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) 
          *dataPtr++ -= *otherDataPtr++;
      } else 
        for (int i=0; i<N; i++) 
          *dataPtr++ += *otherDataPtr++ * otherFact;

    } else if (thisFact == 0.0) {
        // want: this = other * otherFact
        double *dataPtr = values;
        const double *otherDataPtr = other.values;
        if (otherFact == 1.0) {
          // no point doing a multiplication if otherFact == 1.0
          for (int i=0; i<N; i++) 
            *dataPtr++ = *otherDataPtr++;
        } else if (otherFact == -1.0) {
          // no point doing a multiplication if otherFact == 1.0
          for (int i=0; i<N; i++) 
            *dataPtr++ = -(*otherDataPtr++);
        } else 
          for (int i=0; i<N; i++) 
            *dataPtr++ = *otherDataPtr++ * otherFact;
    } else {
      // want: this = this * thisFact + other * otherFact
      double *dataPtr = values;
      const double *otherDataPtr = other.values;
      if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) {
          double value = *dataPtr * thisFact + *otherDataPtr++;
          *dataPtr++ = value;
        }
      } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) {
          double value = *dataPtr * thisFact - *otherDataPtr++;
          *dataPtr++ = value;
        }
      } else 
        for (int i=0; i<N; i++) {
          double value = *dataPtr * thisFact + *otherDataPtr++ * otherFact;
          *dataPtr++ = value;
      }
    }

    // successfull
    return 0;
  }


  template <int NC>
  inline int
  addMatrixVector(double thisFact, const MatrixND<N, NC, double> &m, const Vector &v, double otherFact)
  {
    // check the sizes are compatable
    assert(NC == v.sz);

    // see if quick return
    if (thisFact == 1.0 && otherFact == 0.0)
      return 0;

    else {
      int incr = 1,
             i = N,
             n = NC;
       DGEMV("N", &i, &n,
             &otherFact,
             &m.values[0][0], &i,
             v.theData, &incr,
             &thisFact,
             values,   &incr);
      // successfull
      return 0;
    } 
  }

  template <int NR>
  inline int
  addMatrixTransposeVector(double thisFact, const MatrixND<NR, N, double> &m, const Vector &v, double otherFact)
  {
    // check the sizes are compatable
    assert(NR == v.sz);

    // see if quick return
    if (thisFact == 1.0 && otherFact == 0.0)
      return 0;

    else {
      int incr = 1,
             i = NR,
             n = N;
      DGEMV("T", &i, &n,
            &otherFact,
            &m.values[0][0], &i,
            v.theData, &incr,
            &thisFact,
            values,   &incr);
      // successfull
      return 0;
    } 
  }



  inline int
  addMatrixVector(const double thisFact, const Matrix &m, const Vector &v, const double otherFact)
  {
    // check the sizes are compatable
    assert(N == m.noRows());
    assert(m.noCols() == v.sz);

    // see if quick return
    if (thisFact == 1.0 && otherFact == 0.0)
      return 0;

#ifdef VECTOR_BLAS
    else if (v.sz > 10) {
      int incr = 1,
             i = m.numRows,
             n = m.numCols;
      return
        DGEMV("N", &i, &n,
              &otherFact,
              m.data, &i,
              v.theData, &incr,
              &thisFact,
              values,   &incr);
    }
#endif

    else if (thisFact == 1.0) {

      // want: this += m * v * otherFact
      if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
        int otherSize = v.sz;
        double *matrixDataPtr = m.data;
        const double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          double otherData = *otherDataPtr++;
          for (int j=0; j < N; j++)
            values[j] += *matrixDataPtr++ * otherData;
        }
      } 
      else if (otherFact == -1.0) { // no point doing multiplication if otherFact = -1.0
        int otherSize = v.sz;
        double *matrixDataPtr = m.data;
        const double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          double otherData = *otherDataPtr++;
          for (int j=0; j < N; j++)
            values[j] -= *matrixDataPtr++ * otherData;
        }
      } 
      else { // have to do the multiplication
        int otherSize = v.sz;
        double *matrixDataPtr = m.data;
        double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          double otherData = *otherDataPtr++ * otherFact;
          for (int j=0; j < N; j++)
            values[j] += *matrixDataPtr++ * otherData;
        }
      }
    }

    else if (thisFact == 0.0) {
      
      // want: this = m * v * otherFact
      for (int i=0; i < N; i++)
        values[i] = 0.0;

      if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
        int otherSize = v.sz;
        const double *matrixDataPtr = m.data;
        const double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          const double otherData = *otherDataPtr++;
          for (int j=0; j < N; j++)
            values[j] += *matrixDataPtr++ * otherData;
        }
      } 
      else if (otherFact == -1.0) { // no point doing multiplication if otherFact = -1.0
        int otherSize = v.sz;
        double *matrixDataPtr = m.data;
        const double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          const double otherData = *otherDataPtr++;
          for (int j=0; j < N; j++)
            values[j] -= *matrixDataPtr++ * otherData;
        }
      } else {
        int otherSize = v.sz;
        double *matrixDataPtr = m.data;
        double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          double otherData = *otherDataPtr++ * otherFact;
          for (int j=0; j < N; j++)
            values[j] += *matrixDataPtr++ * otherData;
        }
      }
    }

    else {

      // want: this = this * thisFact + m * v * otherFact
      for (int i=0; i<N; i++)
        values[i] *= thisFact;

      if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
        int otherSize = v.sz;
        double *matrixDataPtr = m.data;
        double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          double otherData = *otherDataPtr++;
          for (int j=0; j < N; j++)
            values[j] += *matrixDataPtr++ * otherData;
        }
      } else if (otherFact == -1.0) { // no point doing multiplication if otherFact = 1.0
        int otherSize = v.sz;
        double *matrixDataPtr = m.data;
        double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          double otherData = *otherDataPtr++;
          for (int j=0; j < N; j++)
            values[j] -= *matrixDataPtr++ * otherData;
        }
      } else {
        int otherSize = v.sz;
        double *matrixDataPtr = m.data;
        double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          double otherData = *otherDataPtr++ * otherFact;
          for (int j=0; j < N; j++)
            values[j] += *matrixDataPtr++ * otherData;
        }
      }
    }
    
    // successfull
    return 0;
  }

  template<typename VecT> inline
  constexpr T
  dot(const VecT &other) const {
    T sum = 0.0;
    for (index_t i = 0; i < N; ++i) {
      sum += values[i] * other[i];
    }
    return sum;
  }

  // Tensor product, also known as the "bun" product
  template <int nc>
  inline OpenSees::MatrixND<N,nc,double>
  bun(const VectorND<nc> &other) const {
    OpenSees::MatrixND<N,nc,double> prod;

    for (int j = 0; j < 3; ++j)
      for (int i = 0; i < 3; ++i)
        prod(i,j) = values[i] * other.values[j];

    return prod;
  }

  constexpr T
  norm() const {
    return sqrt(this->dot(*this));
  }
  
  inline double normalize() {
    double n = norm();

    for (int i=0; i<N; i++)
      values[i] /= n;

    return n;
  }

  template<typename VecT> inline
  VectorND<N> &operator=(const VecT &right) {
    for (int i=0; i< N; i++)
      values[i] = right[i];
    return *this;
  }

  inline
  VectorND<N> &operator/=(const double &right) {
    for (int i=0; i< N; i++)
      values[i] /= right;
    return *this;
  }

  inline
  VectorND<N>  operator/(const double &right) const {
    VectorND<N> res(*this);
    res /= right;
    return res;
  }

  inline
  VectorND<N> &operator*=(const double &right) {
    for (int i=0; i< N; i++)
      values[i] *= right;
    return *this;
  }

  inline
  VectorND<N>  operator*(const double &right) const {
    VectorND<N> res(*this);
    res *= right;
    return res;
  }

  inline
  VectorND<N> &operator+=(const VectorND<N> &right) {
    for (int i=0; i< N; i++)
      values[i] += right[i];
    return *this;
  }

  VectorND<N> &operator+=(const Vector &right) {
    assert(right.Size() == N);
    for (int i=0; i< N; i++)
      values[i] += right[i];
    return *this;
  }

  template <class VecT>
  VectorND<N> operator+(const VecT &right) const {
    VectorND<N> res {*this};
    res += right;
    return res;
  }

  template <class VecT>
  VectorND<N> &operator-=(const VecT &right) {
    for (int i=0; i< N; i++)
      values[i] -= right[i];
    return *this;
  }

  template <class VecT>
  VectorND<N> operator-(const VecT &right) const {
    VectorND<N> res {*this};
    res -= right;
    return res;
  }

  friend std::ostream &
  operator<<(std::ostream &out, const VectorND &vec) {
    out << "{";
    for (int r=0; r<N; r++){
        out << vec[r] << ( r < N-1? ", ": "");
    }
    return out << "}\n";
  }

    /**
    Returns the cross product this vector with another vector.
    @param b the other vector
    @return the cross product.
    */
    template <class VecB, class VecC> inline 
    void cross(const VecB& b, VecC& c) const requires(N==3) {
        c[0] = values[1] * b[2] - values[2] * b[1];
        c[1] = values[2] * b[0] - values[0] * b[2];
        c[2] = values[0] * b[1] - values[1] * b[0];
        return c;
    }


    template <class Vec3T> inline
    VectorND<N> cross(const Vec3T& b) const requires(N==3) {
        VectorND<3> c;
        c[0] = values[1] * b[2] - values[2] * b[1];
        c[1] = values[2] * b[0] - values[0] * b[2];
        c[2] = values[0] * b[1] - values[1] * b[0];
        return c;
    }
};
}

template<int N>
inline OpenSees::VectorND<N>
operator * (double a, const OpenSees::VectorND<N>& b) {
    return b * a;
}
#endif

