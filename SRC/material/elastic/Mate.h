//
//
//
#pragma once
#include <Matrix3D.h>
#include <MatrixSD.h>
#include <MatrixND.h>
#include <TaggedObject.h>
using OpenSees::Matrix3D;
using OpenSees::MatrixND;
using OpenSees::MatrixSD;

class Response;
class Information;

namespace OpenSees {

template <int ndim>
class Mate : public TaggedObject {
  public:

  Mate(int tag) : TaggedObject(tag) {}

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
  virtual const MatrixSD<ndim>& getStress() = 0;

  virtual void setStrain(const MatrixND<ndim,ndim>&) = 0;

#if 0
  virtual MatrixSD<ndim> get_tangent()
  {
  }
#endif
};
}

