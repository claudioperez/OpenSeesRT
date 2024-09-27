#include <Matrix3D.h>
#include <MatrixSD.h>
#include <MatrixND.h>
using OpenSees::Matrix3D;
using OpenSees::MatrixND;
using OpenSees::MatrixSD;

template <int ndim>
class Mate {
  public:

  virtual void set_strain(const MatrixND<ndim,ndim>&)  {
  }

  virtual const MatrixSD<ndim>&   get_stress();
//virtual MatrixND<ndim,ndim> get_stress_tensor();
//virtual MatrixND<ndim,ndim> get_stress_matrix();
};


template <int ndim>
class ElasticMaterial : public Mate<ndim> {
  public:
  ElasticMaterial(double E, double nu) 
    : E(E), nu(nu)
  {
  }


  virtual MatrixND<ndim,ndim> get_stress_matrix() final 
  {
    MatrixND<ndim,ndim> S;
    for (int i=0; i<ndim; i++)
      for (int j=0; j<ndim; j++)
        S(i,j) = Smat(i,j);

    return S;
  }


  virtual const MatrixSD<ndim>& get_stress() final {
    return Smat;
  }

  virtual void set_strain(const MatrixND<ndim,ndim>& F) 
  {

    MatrixSD<ndim> Emat;
    Emat.addMatrixTransposeProduct(0.0, F, F, 0.5);
    Emat.addDiagonal(-0.5);

    if constexpr (ndim == 3) {
      static MatrixND<6,6> ddsdde{0.0};


      double tmp = E / (1.0+nu) / (1.0-2.0*nu);
      ddsdde(0, 0) = (1.0-nu) * tmp;
      ddsdde(1, 1) = ddsdde(1, 1);
      ddsdde(2, 2) = ddsdde(1, 1);
      ddsdde(0, 1) = nu * tmp;
      ddsdde(0, 2) = ddsdde(1, 2);
      ddsdde(1, 2) = ddsdde(1, 2);
      ddsdde(1, 0) = ddsdde(1, 2);
      ddsdde(2, 0) = ddsdde(1, 2);
      ddsdde(2, 1) = ddsdde(1, 2);
      ddsdde(3, 3) = (1.0-2.0*nu) * tmp;
      ddsdde(4, 4) = ddsdde(4, 4);
      ddsdde(5, 5) = ddsdde(4, 4);

      Smat.vector.addMatrixVector(0.0, ddsdde, Emat.vector, 1.0);
    }

  }

private:
  double E, nu;
  MatrixSD<ndim> Smat;
};


