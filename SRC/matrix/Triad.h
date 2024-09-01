//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
//
//
#ifndef Triad_H
#define Triad_H
#include <Matrix.h>
#include <Matrix3D.h>
#include <Vector3D.h>

struct Triad {

  enum Basis {E1=1, E2=2, E3=3};

  Triad()
    : e{{0.0,0.0,0.0}, // e1
        {0.0,0.0,0.0}, // e2
        {0.0,0.0,0.0}} // e3
  {  
  }

  Triad(const Matrix &E)
    : e{{E(0,0),E(1,0),E(2,0)}, // e1
        {E(0,1),E(1,1),E(2,1)}, // e2
        {E(0,2),E(1,2),E(2,2)}} // e3
  {  
  }

  Triad(const OpenSees::Matrix3D &E)
    : e{{E(0,0),E(1,0),E(2,0)}, // e1
        {E(0,1),E(1,1),E(2,1)}, // e2
        {E(0,2),E(1,2),E(2,2)}} // e3
  {  
  }

  constexpr inline
  const Vector3D& operator[](int i) const {
    return e[i-1];
  }

  constexpr inline
  double operator() (int i, int j) const {
    return e[i-1][j];
  }

  const Vector3D e[3];
};
#endif // Triad_H
