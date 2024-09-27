#pragma once
#include <Matrix3D.h>
using OpenSees::Matrix3D;

class Mate {
  public:
  virtual void set_strain(const Matrix3D&);
  virtual Matrix3D get_stress();
};
