
#include <Domain.h>

//
//
//
template<int ndim>
PeriElement<ndim>::PeriElement(int tag, PeriDomain<ndim>& pdomain)
  : Element(tag, 0), edomain(nullptr), pdomain(pdomain)
{
}

template <int ndim>
int
PeriElement<ndim>::getNumExternalNodes() const
{
  return pdomain.pts.size();
}

template <int ndim>
int
PeriElement<ndim>::commitState()
{
  return 0;
}


template <int ndim>
int
PeriElement<ndim>::update()
{
  return 0;
}

template <int ndim>
int
PeriElement<ndim>::revertToStart()
{
  return 0;
}

template <int ndim>
int
PeriElement<ndim>::revertToLastCommit()
{
  return 0;
}

template<int ndim>
const Vector &
PeriElement<ndim>::getResistingForce()
{
  return *pforce;
}

template <int ndim>
const Matrix &
PeriElement<ndim>::getTangentStiff()
{
  const int nn = this->getNumExternalNodes();
  static Matrix K(ndim*nn,ndim*nn);
  return K;
}

template <int ndim>
const Matrix &
PeriElement<ndim>::getInitialStiff()
{
  const int nn = this->getNumExternalNodes();
  static Matrix K(ndim*nn,ndim*nn);
  return K;
}

template <int ndim>
const Matrix &
PeriElement<ndim>::getDamp()
{
  const int nn = this->getNumExternalNodes();
  static Matrix K(ndim*nn,ndim*nn);
  return K;
}

template <int ndim>
const Matrix &
PeriElement<ndim>::getMass()
{
  const int nn = this->getNumExternalNodes();
  static Matrix K(ndim*nn,ndim*nn);
  return K;
}

