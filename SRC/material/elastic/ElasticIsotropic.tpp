

#include <MatrixND.h>
#include <OPS_Stream.h>

using namespace OpenSees;

template<int ndim, PlaneType type>
ElasticIsotropic<ndim,type>::ElasticIsotropic(int tag, double E, double nu, double rho)
  : Mate<ndim>(tag), E(E), nu(nu), rho(rho)
{

}


template <int ndim, PlaneType type>
const char*
ElasticIsotropic<ndim,type>::getClassType() const
{
  return "ElasticIsotropic";
}


template <int ndim, PlaneType type>
Mate<ndim>*
ElasticIsotropic<ndim,type>::getCopy()
{
  return new ElasticIsotropic<ndim,type>(this->getTag(), E, nu, rho);
}

template<int ndim, PlaneType type>
const MatrixSD<ndim>& 
ElasticIsotropic<ndim,type>::getStress()
{
  return Smat;
}

template<int ndim, PlaneType type>
int
ElasticIsotropic<ndim,type>::setTrialStrain(const MatrixND<ndim,ndim>& F) 
{
  // Convert deformation gradient F to small strain tensor e
  MatrixSD<ndim,true> e;
  e.addMatrixTransposeProduct(0.0, F, F, 0.5);
  e.addDiagonal(-0.5);
  return setTrialStrain(e);
}

template<int ndim, PlaneType type>
int
ElasticIsotropic<ndim,type>::setTrialStrain(const MatrixSD<ndim,true>& e)
{
  if constexpr (ndim == 3) {
    static MatrixND<6,6> ddsdde{{{0.0}}};

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
    ddsdde(3, 3) = (1.0 - 2.0*nu) * tmp;
    ddsdde(4, 4) = ddsdde(4, 4);
    ddsdde(5, 5) = ddsdde(4, 4);

    Smat.vector = ddsdde * e.vector;
  }
}

template<int ndim, PlaneType type>
MatrixSD<ne>
ElasticIsotropic<ndim, type>::getTangent()
{
  if constexpr (ndim == 3) {
    MatrixND<6,6> ddsdde{{{0.0}}};

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
    ddsdde(3, 3) = (1.0 - 2.0*nu) * tmp;
    ddsdde(4, 4) = ddsdde(4, 4);
    ddsdde(5, 5) = ddsdde(4, 4);

    return ddsdde;
  }
}

template<int ndim, PlaneType type>
void
ElasticIsotropic<ndim,type>::Print(OPS_Stream& s, int flag) 
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "Elastic Isotropic Material Model" << "\n";
        s << "\tE:  "   << E   << "\n";
        s << "\tv:  "   << nu  << "\n";
        s << "\trho:  " << rho << "\n";

    } else if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << OPS_PRINT_JSON_ELEM_INDENT << "{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"" << this->getClassType() << "\", ";
        s << "\"E\": "   << E   << ", ";
        s << "\"nu\": "  << nu  << ", ";
        s << "\"rho\": " << rho;
	s << "}";
    }
}
