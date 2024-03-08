//
// Claudio Perez
//
#ifndef Matrix3D_H
#define Matrix3D_H

#include "MatrixND.h"
#include "Vector3D.h"
#include "routines/SY3.h"

namespace OpenSees {
class  Matrix3D: public MatrixND<3,3> {

public:
  template<class VecT> void addSpin(const VecT& V);
  template<class VecT> void addSpin(const VecT& V, const double scale);
  template<class VecT> void addSpinSquare(const VecT& V, const double scale);
  // inline void addEye(Vector3D<double>& v, double scale);
  // void addDev;

  int symeig(Vector3D<double>& vals) {
    double work[3][3];
    cmx_eigSY3(values, work, vals.values);
    return 0;
  }
};

template< class VecT> inline void 
Matrix3D::addSpin(const VecT& v, double mult)
{
   const double v0 = mult*v[0],
                v1 = mult*v[1],
                v2 = mult*v[2];

  (*this)(0, 0) += 0.0;   (*this)(0, 1) += -v2;     (*this)(0, 2) +=  v1;
  (*this)(1, 0) +=  v2;   (*this)(1, 1) += 0.00;    (*this)(1, 2) += -v0;
  (*this)(2, 0) += -v1;   (*this)(2, 1) +=  v0;     (*this)(2, 2) += 0.0;
}

template <class VecT> inline void
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
}

} // namespace OpenSees


#endif // Matrix3D_H
