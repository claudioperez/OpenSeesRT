//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Claudio Perez
//
#ifndef Matrix3D_H
#define Matrix3D_H

#include "MatrixND.h"
#include "Vector3D.h"
#include "routines/SY3.h"

namespace OpenSees {


using Matrix3D = MatrixND<3,3,double>;

#if 0
class  Matrix3D: public MatrixND<3,3,double> {

public:

  template<class VecT> Matrix3D& addSpin(const VecT& V);
  template<class VecT> Matrix3D& addSpin(const VecT& V, double scale);
  template<class VecT> Matrix3D& addSpinSquare(const VecT& V, const double scale);
  template<class VecT> void addSpinProduct(const VecT& a, const Vector3D& b, const double scale);
  template<class VecT> void addMatrixSpinProduct(const Matrix3D& A, const VecT& b, const double scale);
  template<class MatT> void addSpinMatrixProduct(const Vector3D& a, const MatT& B, const double scale);


  Vector3D operator*(const Vector3D&v);

  Matrix3D operator+(const Matrix3D &B) const {
    const Matrix3D &A = *this;
    return Matrix3D {{{
      { A(0,0)+B(0,0), A(1,0)+B(1,0), A(2,0)+B(2,0) },
      { A(0,1)+B(0,1), A(1,1)+B(1,1), A(2,1)+B(2,1) },
      { A(0,2)+B(0,2), A(1,2)+B(1,2), A(2,2)+B(2,2) }
    }}};
  }

};
#endif

#if 0
inline Vector3D
Matrix3D::operator*(const Vector3D&v)
{
  return Vector3D {{
    (*this)(0,0)*v[0] + (*this)(0,1)*v[1] + (*this)(0,2)*v[2],
    (*this)(1,0)*v[0] + (*this)(1,1)*v[1] + (*this)(1,2)*v[2],
    (*this)(2,0)*v[0] + (*this)(2,1)*v[1] + (*this)(2,2)*v[2]
  }};
}


template< class VecT> inline 
Matrix3D& 
Matrix3D::addSpin(const VecT& v)
{
   const double v0 = v[0],
                v1 = v[1],
                v2 = v[2];

  (*this)(0, 0) += 0.0;   (*this)(0, 1) += -v2;     (*this)(0, 2) +=  v1;
  (*this)(1, 0) +=  v2;   (*this)(1, 1) +=  0.0;    (*this)(1, 2) += -v0;
  (*this)(2, 0) += -v1;   (*this)(2, 1) +=  v0;     (*this)(2, 2) += 0.0;

  return *this;
}


template<class VecT> inline 
Matrix3D&
Matrix3D::addSpin(const VecT& v, double mult)
{
   const double v0 = mult*v[0],
                v1 = mult*v[1],
                v2 = mult*v[2];

  (*this)(0, 0) += 0.0;   (*this)(0, 1) += -v2;     (*this)(0, 2) +=  v1;
  (*this)(1, 0) +=  v2;   (*this)(1, 1) += 0.00;    (*this)(1, 2) += -v0;
  (*this)(2, 0) += -v1;   (*this)(2, 1) +=  v0;     (*this)(2, 2) += 0.0;
  return *this;
}


template <class VecT> inline
Matrix3D& 
Matrix3D::addSpinSquare(const VecT& v, const double scale)
{
  const double v1 = v[0],
               v2 = v[1],
               v3 = v[2];

  (*this)(0,0) += scale*( -v2*v2 - v3*v3 );
  (*this)(1,1) += scale*( -v1*v1 - v3*v3 );
  (*this)(2,2) += scale*( -v1*v1 - v2*v2 );

  (*this)(0,1) += scale*(  v1*v2 );
  (*this)(1,0) += scale*(  v1*v2 );
  (*this)(2,0) += scale*(  v1*v3 );
  (*this)(0,2) += scale*(  v1*v3 );
  (*this)(1,2) += scale*(  v2*v3 );
  (*this)(2,1) += scale*(  v2*v3 );
  return *this;
}


template<class VecT> inline
void Matrix3D::addSpinProduct(const VecT& a, const Vector3D& b, const double scale)
{
  // a^b^ = boa - a.b 1
  // where 'o' denotes the tensor product and '.' the dot product
  //
  this->addTensorProduct(b, a, scale);
  this->addDiagonal(-b.dot(a)*scale);
}

template<class VecT> inline
void Matrix3D::addMatrixSpinProduct(const Matrix3D& A, const VecT& b, const double scale)
{
  // this += s*A*[b^]
  // where b^ is the skew-symmetric representation of the three-vector b, s is a scalar,
  // and A a 3x3 matrix.
  //
  (*this)(0, 0) += scale*( A(0,1)*b[2] - A(0,2)*b[1]);
  (*this)(0, 1) += scale*(-A(0,0)*b[2] + A(0,2)*b[0]);
  (*this)(0, 2) += scale*( A(0,0)*b[1] - A(0,1)*b[0]);
  (*this)(1, 0) += scale*( A(1,1)*b[2] - A(1,2)*b[1]);
  (*this)(1, 1) += scale*(-A(1,0)*b[2] + A(1,2)*b[0]);
  (*this)(1, 2) += scale*( A(1,0)*b[1] - A(1,1)*b[0]);
  (*this)(2, 0) += scale*( A(2,1)*b[2] - A(2,2)*b[1]);
  (*this)(2, 1) += scale*(-A(2,0)*b[2] + A(2,2)*b[0]);
  (*this)(2, 2) += scale*( A(2,0)*b[1] - A(2,1)*b[0]);
}

template<class MatT> inline
void Matrix3D::addSpinMatrixProduct(const Vector3D& a, const MatT& B, const double scale)
{
  // this += s*[a^]*B
  // where a^ is the skew-symmetric representation of the three-vector a, s is a scalar,
  // and B a 3x3 matrix.
  //
  (*this)(0, 0) += scale*( -B(1,0)*a[2] + B(2,0)*a[1]);
  (*this)(0, 1) += scale*( -B(1,1)*a[2] + B(2,1)*a[1]);
  (*this)(0, 2) += scale*( -B(1,2)*a[2] + B(2,2)*a[1]);
  (*this)(1, 0) += scale*(  B(0,0)*a[2] - B(2,0)*a[0]);
  (*this)(1, 1) += scale*(  B(0,1)*a[2] - B(2,1)*a[0]);
  (*this)(1, 2) += scale*(  B(0,2)*a[2] - B(2,2)*a[0]);
  (*this)(2, 0) += scale*( -B(0,0)*a[1] + B(1,0)*a[0]);
  (*this)(2, 1) += scale*( -B(0,1)*a[1] + B(1,1)*a[0]);
  (*this)(2, 2) += scale*( -B(0,2)*a[1] + B(1,2)*a[0]);
}
#endif

} // namespace OpenSees


#endif // Matrix3D_H
