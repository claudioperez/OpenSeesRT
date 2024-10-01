#include "Mate.h"

template <int ndim>
class ElasticIsotropic : public Mate<ndim> {
  public:
  ElasticIsotropic(double E, double nu);

  virtual Mate<ndim>* getCopy() final;
//virtual MatrixND<ndim,ndim> get_stress_matrix() final;
  virtual const MatrixSD<ndim>& get_stress() final;
  virtual void set_strain(const MatrixND<ndim,ndim>& F) final;

private:
  double E, nu;
  MatrixSD<ndim> Smat;
};

#include "ElasticIsotropic.tpp"
