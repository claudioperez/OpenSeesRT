#include <Mate.h>
#include <NosbBase.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <Vector3D.h>
#include <Matrix3D.h>

using namespace OpenSees;

template <int ndim,int mfam>
class NosbProj : public NosbBase<ndim> {
public :
  int numfam;
  double scale;
  Matrix3D Kmat;
  PeriParticle<ndim>& center;
  std::array<PeriParticle<ndim>*, mfam> family;
  std::array<Mate*,               mfam> materials;
  std::array<double,              mfam> omega;

public:
  NosbProj(PeriParticle<ndim>* center, PeriDomain<ndim>& domain, Mate* material)
    : center(*center), numfam(center->numfam) 
  {
    for (int i = 0; i<numfam; i++) {
      family[i] = &domain.pts[center->nodefam[i]];
      materials[i] =  material;
    }
  }


  double
  bond_omega(int i)
  {
//  const PeriParticle<3>& other = 
    // TODO
//  omega[i] = 1.0;
    return 1.0; // omega[i];
  }

  void
  init_shape()
  {

    // First populate omega
    for (int i=0; i<numfam; i++) {
      omega[i] = bond_omega(i)/center.vol_h;
    }

    // Integrate K
    Kmat.zero();

    for (int i = 0; i<numfam; i++) {
      const Vector3D xi = family[i]->coord - center.coord;
      Kmat +=  xi.bun(xi) * center.vol[i] * omega[i];
    }

    if (ndim == 2) 
      Kmat(2, 2) = 1.0;

    // Invert K
    Kmat.invert();
  }

  MatrixND<ndim,ndim>
  getA(const VectorND<ndim>& xi)
  {
    // Initialize with zeros
    MatrixND<ndim,ndim> A{0.0};

    A.addDiagonal(1.0);

    A.addTensorProduct(xi, xi, -1.0/xi.dot(xi));

    return A;
  }


  VectorND<ndim>
  get_T2(int i, const VectorND<ndim>& xi)
  {
    const MatrixND<ndim,ndim>&   P = materials[i]->get_stress();
    VectorND<ndim>               B = P * xi;
    B /= xi.dot(xi);
    return B;
  }


  //
  Matrix3D
  getKinv()
  {
    Matrix3D Kinv;
    // TODO
    return Kinv;
  }


  void
  form_trial()
  {
    // update materials (materials)
    const Matrix3D Kinv = this->getKinv();

    Matrix3D Nmat{0.0};
    for (int i = 0; i < numfam; i++) {
      PeriParticle<ndim>& other = *family[i];
      const Vector3D eta  = other.disp  - center.disp;
      const Vector3D xi   = other.coord - center.coord;
      const Vector3D zeta = xi + eta;
      Nmat += omega[i] * zeta.bun(xi) * center.vol[i];
    }

    Matrix3D Fmat = Nmat * Kinv;

    if (ndim == 2)
      Fmat(2, 2) = 1.0;

    for (int i = 0; i < numfam; i++) {
      PeriParticle<ndim>& other = *family[i];
      const Vector3D eta  = other.disp - center.disp;
      const Vector3D xi   = other.coord - center.coord;
      const Vector3D zeta = xi + eta;

      const Matrix3D& Amat = this->getA(xi);
      const Matrix3D  Kmat = zeta.bun(xi) / xi.dot(xi);
      // Fmat <- Fmat * (I - xi\otimes xi / |xi|^2) + zeta\otimes xi / |xi|^2
      materials[i]->set_strain(Fmat*Amat + Kmat);
    }

  }

public:
  virtual MatrixND<ndim,ndim>
  sum_PKinv()
  {

    const Matrix3D Kinv = this->getKinv();
    Matrix3D Qmat{0.0};
    for (int i = 0; i < numfam; i++) {
      PeriParticle<ndim>& other = *family[i];
      const Vector3D xi = other.coord - center.coord;
      const Matrix3D& P = materials[i]->get_stress();
      const Matrix3D& A = this->getA(xi);
      Qmat.addMatrix(P*A*Kinv, omega[i]*center.vol[i]);
    }
    return Qmat;
  }

  virtual VectorND<ndim>
  bond_force(int i, MatrixND<ndim,ndim>& Qmat)
  {
      VectorND<ndim> pforce;
      const Matrix3D Kinv = this->getKinv();

      PeriParticle<ndim>& other = *family[i];

      const Vector3D xi = other.coord - center.coord;
      Vector3D B = this->get_T2(i, xi);

      Vector3D  T_j = B + Qmat*xi;

      double w = omega[i] * center.vol[i] * center.correction[i];

      return T_j*w;
  }
};


