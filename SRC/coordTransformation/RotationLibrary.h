//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
//
#pragma once
#include <VectorND.h>
#include <Vector3D.h>
#include <MatrixND.h>
#include <Matrix3D.h>
using OpenSees::Matrix3D;
using OpenSees::MatrixND;
using OpenSees::VectorND;

template <class VecT> static inline void
getTangScaledPseudoVectorFromQuaternion(const VectorND<4> &q, VecT& w)
{
  const double q0 = q[3];

  for (int i = 0; i < 3; i++)
    w[i] = 2.0 * q[i]/q0;

  return ;
}

static inline VectorND<4>
quaternionProduct(const VectorND<4> &qa, const VectorND<4> &qb)
{
    const double qa0 = qa[0],
                 qa1 = qa[1],
                 qa2 = qa[2],
                 qa3 = qa[3],
                 qb0 = qb[0],
                 qb1 = qb[1],
                 qb2 = qb[2],
                 qb3 = qb[3];

    // calculate the dot product qa.qb
    const double qaTqb = qa0*qb0 + qa1*qb1 + qa2*qb2;

    // calculate the cross-product qa x qb
    const double
      qaxqb0 = qa1*qb2 - qa2*qb1,
      qaxqb1 = qa2*qb0 - qa0*qb2,
      qaxqb2 = qa0*qb1 - qa1*qb0;

    // calculate the quaternion product
    VectorND<4> q12;
    q12[0] = qa3*qb0 + qb3*qa0 - qaxqb0;
    q12[1] = qa3*qb1 + qb3*qa1 - qaxqb1;
    q12[2] = qa3*qb2 + qb3*qa2 - qaxqb2;
    q12[3] = qa3*qb3 - qaTqb;
    return q12;
}

static inline Vector3D
LogC90(const MatrixND<3,3> &R)
{
  return Vector3D {
    std::asin(0.5*(R(1,2) - R(2,1))),
    std::asin(0.5*(R(0,1) - R(1,0))),
    std::asin(0.5*(R(0,2) - R(2,0))),
  };
}

static inline Matrix3D
CaySO3(const Vector3D &w)
{
  // Cayley map: for a rotation matrix given the tangent-scaled pseudo-vector

  // R = I + (S + S*S/2)/(1 + w' * w / 4);
  const double c = 1.0/(1 + w.dot(w)/4.0);

  Matrix3D R;
  R.zero();
  R.addDiagonal(1.0);
  R.addSpin(w, c);
  R.addSpinSquare(w, 0.5*c);

  return R;
}

static inline void 
getRotationMatrixFromQuaternion(const VectorND<4> &q, Matrix3D& R)
{
    // R = (q0^2 - q' * q) * I + 2 * q * q' + 2*q0*S(q);

    const double factor = q[3]*q[3] - (q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);

    R.zero();

    for (int i = 0; i < 3; i++)
      R(i,i) = factor;

    R.addTensorProduct(q, q, 2.0);
    R.addSpin(q, 2.0*q[3]);
}


static inline VectorND<4>
VersorFromMatrix(const Matrix3D &R)
{
    VectorND<4> q;
    // obtains the normalised quaternion from the rotation matrix
    // using Spurrier's algorithm

    const double trR = R(0,0) + R(1,1) + R(2,2);

    // a = max ([trR R(0,0) R(1,1) R(2,2)]);
    double a = trR;
    for (int i = 0; i < 3; i++)
      if (R(i,i) > a)
        a = R(i,i);

    if (a == trR) {
      q[3] = sqrt(1+a)*0.5;

      for (int i = 0; i < 3; i++) {
        int j = (i+1)%3;
        int k = (i+2)%3;
        q[i] = (R(k,j) - R(j,k))/(4*q[3]);
      }
    }
    else {
      for (int i = 0; i < 3; i++)
        if (a == R(i,i)) {
          int j = (i+1)%3;
          int k = (i+2)%3;

          q[i] = sqrt(a*0.5 + (1 - trR)/4.0);
          q[3] = (R(k,j) - R(j,k))/(4*q[i]);
          q[j] = (R(j,i) + R(i,j))/(4*q[i]);
          q[k] = (R(k,i) + R(i,k))/(4*q[i]);
        }
    }
    return q;
}

static inline VectorND<4>
getQuaternionFromPseudoRotVector(const Vector  &theta)
{
    VectorND<4> q;      // normalized quaternion

    double t = theta.Norm();
    if (t == 0)
        q.zero();

    else {
        const double factor = sin(t*0.5)/ t;
        for (int i = 0; i < 3; i++)
            q[i] = theta(i) * factor;
    }

    q[3] = cos(t*0.5);

    return q;
}

