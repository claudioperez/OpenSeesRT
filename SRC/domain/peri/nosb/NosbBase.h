#pragma once
#include <MatrixND.h>
#include <VectorND.h>
using namespace OpenSees;

template <int ndim>
class NosbBase {
  public:
  virtual void           init_shape();
  virtual void           form_trial();

  virtual MatrixND<ndim,ndim> sum_PKinv();
  virtual VectorND<ndim> bond_force(int i, MatrixND<ndim,ndim>& Q);
};

