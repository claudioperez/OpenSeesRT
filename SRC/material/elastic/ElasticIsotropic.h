#include "Mate.h"

enum class PlaneType {None, Stress, Strain};

template <int ndim, PlaneType = PlaneType::None>
class ElasticIsotropic : public Mate<ndim> {
  public:
  ElasticIsotropic(int tag, double E, double nu, double rho);

  virtual Mate<ndim>* getCopy() final;
  virtual const char* getClassType() const final;
//virtual MatrixND<ndim,ndim> get_stress_matrix() final;
  virtual const MatrixSD<ndim>& get_stress() final;
  virtual void set_strain(const MatrixND<ndim,ndim>& F) final;

  virtual void Print(OPS_Stream& s, int flag);

private:
  double E, nu, rho;
  MatrixSD<ndim> Smat;
};

#include "ElasticIsotropic.tpp"
