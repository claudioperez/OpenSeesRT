
template <int ndim>
ElasticIsotropic<ndim>::ElasticIsotropic(double E, double nu)
	: E(E), nu(nu)
{
}

template <int ndim>
Mate<ndim> *
ElasticIsotropic<ndim>::getCopy()
{
	return new ElasticIsotropic<ndim>(E, nu);
}

template <int ndim>
const MatrixSD<ndim> &
ElasticIsotropic<ndim>::get_stress()
{
	return Smat;
}

template <int ndim>
void ElasticIsotropic<ndim>::set_strain(const MatrixND<ndim, ndim> &F)
{
	MatrixSD<ndim> Emat{0.0};
	if (ndim == 2){
		Emat(0, 0) = F(0, 0) - 1.0;
		Emat(1, 1) = F(1, 1) - 1.0;
		Emat(0, 1) = 0.5 * F(0, 1) + 0.5 * F(1, 0);
	}
	// Emat.addMatrixTransposeProduct(0.0, F, F, 0.5);
	// Emat.addDiagonal(-0.5);

	if constexpr (ndim == 3)
	{
		static MatrixND<6, 6> ddsdde{0.0};

		double tmp = E / (1.0 + nu) / (1.0 - 2.0 * nu);
		ddsdde(0, 0) = (1.0 - nu) * tmp;
		ddsdde(1, 1) = ddsdde(1, 1);
		ddsdde(2, 2) = ddsdde(1, 1);
		ddsdde(0, 1) = nu * tmp;
		ddsdde(0, 2) = ddsdde(1, 2);
		ddsdde(1, 2) = ddsdde(1, 2);
		ddsdde(1, 0) = ddsdde(1, 2);
		ddsdde(2, 0) = ddsdde(1, 2);
		ddsdde(2, 1) = ddsdde(1, 2);
		ddsdde(3, 3) = (1.0 - 2.0 * nu) * tmp;
		ddsdde(4, 4) = ddsdde(4, 4);
		ddsdde(5, 5) = ddsdde(4, 4);

		Smat.vector.addMatrixVector(0.0, ddsdde, Emat.vector, 1.0);
	}
	else if constexpr (ndim == 2)
	{
		static MatrixND<3, 3> ddsdde{0.0};
		// this corresponds to plane strain case
		const double tmp = E / (1.0 + nu) / (1.0 - 2.0 * nu);
		ddsdde(0, 0) = (1.0 - nu) * tmp;
		ddsdde(1, 1) = (1.0 - nu) * tmp;
		ddsdde(0, 1) = nu * tmp;
		ddsdde(1, 0) = nu * tmp;
		ddsdde(2, 2) = (1.0 - 2.0 * nu) * tmp;
		const VectorND<3> Evec = Emat.vector; // I have no idea why this is necessary
		Smat.vector.addMatrixVector(0.0, ddsdde, Emat.vector, 1.0);
	}
}
