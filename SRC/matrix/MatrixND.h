/* 
 * Objectives:
 * - little to no overhead above C-style arrays
 * - value semantics; objects do not decay to pointers;
 *
 * This code is influenced by the following sources
 *  list initialization:
 *  - https://stackoverflow.com/questions/42068882/list-initialization-for-a-matrix-class
 *
 *  style/practices
 *  - https://quuxplusone.github.io/blog/2021/04/03/static-constexpr-whittling-knife/
 * 
 *  Operator overloading / semantics
 *  - https://stackoverflow.com/questions/9851188/does-it-make-sense-to-use-move-semantics-for-operator-and-or-operator/9851423#9851423
 *
 *  compile-time template restrictions/concepts:
 *  - https://codereview.stackexchange.com/questions/259038/compile-time-matrix-class
 *    (C++ 20)
 *  - https://github.com/calebzulawski/cotila/
 *    (C++ 17)
 */

// Claudio Perez
#ifndef MatrixND_H
#define MatrixND_H
#include <math.h>
#include <assert.h>
#include <array>
#include <stdexcept>
#include <functional>
#include <iostream> // overloading <<


#if __cplusplus < 202000L
#define consteval
#define requires(X)
#endif

#define G23_STACK_MAX 10

namespace OpenSees {

template <index_t NR, index_t NC, typename T=double>
requires(NR > 0 && NC > 0)
struct MatrixND {
  double values[NC][NR];

  constexpr std::array<T, NC> &
  operator[](index_t index) {return values[index];}

  constexpr const std::array<T, NC> & // [i] indexing
  operator[](index_t index) const {return values[index];}
  
  constexpr T & // (i,j) indexing
  operator()(index_t index_r, index_t index_c) {
#ifdef BOUNDS_CHECK
  assert(index_r >= 0 && index_c >= 0);
  assert(index_r < NR && index_c < NC);
#endif
    // column-major
    // return values[index_c*NR + index_r];
    return values[index_c][index_r];
  }

  constexpr const T & // (i,j) indexing
  operator()(index_t index_r, index_t index_c) const {
    if (index_r < 0 || index_c < 0) {
      throw std::out_of_range("negative MatrixND index");
    } else if (index_r >= NR || index_c >= NC) {
      throw std::out_of_range("compound access out of range");
    }
    return values[index_c][index_r];
  }

  constexpr VectorND<NR>
  column(index_t index) const {
    assert(index >= 0);
    assert(index < NC);

    // TODO: ugly temporary implementation
    return *(VectorND<NR>*)(&(values[index]));
  }

  constexpr VectorND<NC>
  row(index_t index) const {
    assert(index >= 0);
    assert(index < NR);

    VectorND<NC,T> rw;
    for (index_t j = 0; j < NC; ++j) {
      rw[j] = values[j][index];
    }
    return rw;
  }

  consteval VectorND<2,int>
  size() const {return {NR, NC};}
  
  constexpr MatrixND &
  operator+=(const double value) {
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        values[j][i] += value;
      }
    }
    return *this;
  }

  constexpr MatrixND &
  operator+=(const MatrixND &other) {
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        values[j][i] += other.values[j][i];
      }
    }
    return *this;
  }
  
  constexpr MatrixND &
  operator-=(const MatrixND &other) {
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        values[j][i] -= other.values[j][i];
      }
    }
    return *this;
  }

  constexpr MatrixND &
  operator*=(T const scalar) {
    for (index_t j = 0; j < NC; ++j){
      for (index_t i = 0; i < NR; ++i) {
        values[j][i] *= scalar;
      }
    }
    return *this;
  }

  constexpr MatrixND &
  operator/=(T const scalar) {
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        values[j][i] /= scalar;
      }
    }
    return *this;
  }


// Notes on operators:
// - define friend operators inside class for use as header-only library
// - LHS is passed by copied value
  friend constexpr MatrixND
  operator+(MatrixND left, const MatrixND &right) {left += right; return left;}

  friend constexpr MatrixND
  operator-(MatrixND left, const MatrixND &right) {left -= right; return left;}
  
  friend constexpr MatrixND // scalar * Matrix
  operator*(T scalar, MatrixND mat) {mat *= scalar; return mat;}

  friend constexpr MatrixND // Matrix * Matrix
  operator*(MatrixND mat, T scalar) {mat *= scalar; return mat;}
  
  template <index_t J>
  constexpr friend MatrixND<NR, J>
  operator*(const MatrixND<NR, NC> &left, const MatrixND<NC, J> &right) {
    MatrixND<NR, J> prod;
    for (index_t i = 0; i < NR; ++i) {
      for (index_t j = 0; j < J; ++j) {
        for (index_t k = 0; k < NC; ++k) {
          prod(i, j) += left(i,k) * right(k,j);
        }
      }
    }
    return prod;
  }

  friend constexpr MatrixND
  operator/(MatrixND mat, T scalar) {mat /= scalar; return mat;}
  
  friend std::ostream &
  operator<<(std::ostream &out, MatrixND const &mat) {
    out << "{";
    for (int r=0; r<NR; r++){
      out << "{";
      for (int c=0; c<NC; c++)
        out << mat(r,c);
      out << "}, ";
    }
    return out << "}\n";
  }

/**********************
 *
 */
  constexpr MatrixND<NC, NR>
  transpose() const {
    MatrixND<NC, NR> result = {};
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        result.values[i][j] = values[j][i];
      }
    }
    return result;
  }

  constexpr T
  trace() const
  requires(NR == NC) 
  {
    T sum = 0.0;
    for (index_t i = 0; i < NR; ++i) {
      sum += values[i][i];
    }
    return sum;
  }

  constexpr T
  determinant() const
  requires(NR == NC && NR == 2)
  {
    return values[0][0] * values[1][1] - values[0][1] * values[1][0];
  }

};


template <int nr, int nc, typename T=double>
MatrixND(const T (&)[nc][nr])->MatrixND<nr, nc, T>;

namespace linalg {

template<int nr, int nc, typename T=double>
inline MatrixND<nr, nc, T> zeros(void){
  MatrixND<nr, nc, T> m;
  for (int i=0; i<nr; i++)
    for (int j=0; j<nc; j++)
      m(i,j) = T(0.0);
  return m;
}


template<int nr, int nc>
// constrain nr and nc to avoid blowing out the stack
requires (nr < G23_STACK_MAX && nc < G23_STACK_MAX)
inline MatrixND<nr,nc>
solve(MatrixND<nr,nr> A, MatrixND<nr,nc> B){
  int pivot_ind[nr] ;
  MatrixND<nr,nc> X = B; // X will be overwritten with the solution
  // LAPACKE_dgesv(LAPACK_COL_MAJOR, nr, nc, &A(0,0), nr, pivot_ind, &X(0,0), nc);
  return X;
}

} // namespace linalg

} // namespace g23

#endif // MatrixND_H

