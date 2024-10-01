#pragma once
#include <MatrixND.h>
#include <VectorND.h>
using namespace OpenSees;

template <int ndim>
class NosbBase
{
public:
  virtual ~NosbBase() {};
	// ============================================
	// MEMBER FUNCTIONS
	// ============================================
	virtual void init_shape() = 0;
	virtual void form_trial() = 0;
	virtual VectorND<ndim> bond_force(int i, MatrixND<ndim, ndim> &Q) = 0;
	virtual MatrixND<ndim, ndim> sum_PKinv() = 0;

	// ============================================
	// MEMBER DATA
	// ============================================
	
};
