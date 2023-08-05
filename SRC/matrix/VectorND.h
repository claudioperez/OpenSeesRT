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
#ifndef VectorND_H
#define VectorND_H
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

namespace OpenSeees {
typedef int index_t;

template <index_t N, typename T=double> 
requires(N > 0)
struct VectorND {
  T values[N];

  constexpr T&
  operator[](index_t index) {return values[index];}
  constexpr const T&
  operator[](index_t index) const {return values[index];}

  constexpr T
  dot(const VectorND<N> &other) const {
    T sum = 0.0;
    for (index_t i = 0; i < N; ++i) {
      sum += values[i] * other.values[i];
    }
    return sum;
  }

  friend std::ostream &
  operator<<(std::ostream &out, VectorND const &vec) {
    out << "{";
    for (int r=0; r<N; r++){
        out << vec[r] << ( r < N-1? ", ": "");
    }
    return out << "}\n";
  }
};
}
#endif
