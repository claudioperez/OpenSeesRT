#ifndef Matrix3D_H
#define Matrix3D_H
#include "MatrixND.h"
#include "Vector3D.h"
#include "routines/SY3.h"

namespace OpenSees {
class  Matrix3D: public MatrixND<3,3> {
  public:
  int symeig(Vector3D<double>& vals) {
    double work[3][3];
    cmx_eigSY3(values, work, vals.mData);
    return 0;
  };
};
}
#endif // Matrix3D_H
