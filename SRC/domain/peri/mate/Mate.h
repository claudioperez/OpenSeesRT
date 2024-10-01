//
//
//
#pragma once
#include <Matrix3D.h>
#include <MatrixSD.h>
#include <MatrixND.h>
using OpenSees::Matrix3D;
using OpenSees::MatrixND;
using OpenSees::MatrixSD;

class Response;
class Information;

template <int ndim>
class Mate {
  public:

  virtual Mate<ndim>* getCopy() = 0;

  virtual int commit()
  {
    return 0;
  }

  virtual int revertToStart()
  {
    return 0;
  }

  virtual int revertToLast()
  {
    return 0;
  }
#if 0
  virtual Response *setResponse (const char **argv, int argc, OPS_Stream &s);
  virtual int getResponse (int responseID, Information &matInformation);
#endif
  virtual const MatrixSD<ndim>& get_stress() = 0;

  virtual void set_strain(const MatrixND<ndim,ndim>&) = 0;

#if 0
  virtual MatrixSD<ndim> get_tangent()
  {
  }
#endif
};

