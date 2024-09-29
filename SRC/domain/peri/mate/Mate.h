#pragma once
#include <Matrix3D.h>
#include <MatrixSD.h>
#include <MatrixND.h>
using OpenSees::Matrix3D;
using OpenSees::MatrixND;
using OpenSees::MatrixSD;

template <int ndim>
class Mate {
  public:

  virtual const MatrixSD<ndim>& get_stress() = 0;

  virtual void set_strain(const MatrixND<ndim,ndim>&) = 0;

};

