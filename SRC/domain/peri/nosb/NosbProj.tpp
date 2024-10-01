#include <Mate.h>
#include <NosbBase.h>
#include <Logging.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <Vector3D.h>
#include <Matrix3D.h>


template <int ndim, int maxfam>
NosbProj<ndim, maxfam>::NosbProj(PeriParticle<ndim> *center, PeriDomain<ndim> &domain, Mate<ndim> *material)
    : center(center), numfam(center->numfam)
{
    for (int i = 0; i < numfam; i++)
    {
        neigh[i] = &domain.pts[center->nodefam[i]];
        materials[i] = material->getCopy();
    }
}

template <int ndim, int maxfam>
double 
NosbProj<ndim, maxfam>::bond_omega(int i)
{
    // in this formulation, omega is constant 
	// for computing all the bond-associated quantities
    
    return 1.0; // omega[i];
}

template <int ndim, int maxfam>
void 
NosbProj<ndim, maxfam>::init_shape()
{
    // First populate omega
    for (int i = 0; i < numfam; i++)
    {
        omega[i] = bond_omega(i) / center->vol_h;
    }

    // Integrate K
    // note that Kinv is invariant in this formulation, but may not be in general
    Kinv.zero();

    for (int i = 0; i < numfam; i++)
    {
        const VectorND<ndim> xi = neigh[i]->coord - center->coord;
        Kinv += omega[i] * xi.bun(xi) * center->vol[i];
    }
    opserr << Matrix(Kinv) << "\n";

    // Invert K
    Kinv.invert();
    opserr << Matrix(Kinv) << "\n";
}

template <int ndim, int maxfam>
MatrixND<ndim, ndim>
NosbProj<ndim, maxfam>::get_A(const VectorND<ndim> &xi)
{
    // Initialize with zeros
    MatrixND<ndim, ndim> A{{0.0}};

    A.addDiagonal(1.0);

    A.addTensorProduct(xi, xi, -1.0 / xi.dot(xi));

    return A;
}

template <int ndim, int maxfam>
VectorND<ndim>
NosbProj<ndim, maxfam>::get_T2(int i, const VectorND<ndim> &xi)
{
    const MatrixSD<ndim>& P = materials[i]->get_stress();
    VectorND<ndim> B = P * xi;
    B /= xi.dot(xi);
    return B;
}

template <int ndim, int maxfam>
void
NosbProj<ndim, maxfam>::form_trial()
{
    // update materials (materials)

    MatrixND<ndim, ndim> Nmat{{0.0}};

    for (int i = 0; i < numfam; i++)
    {
        // PeriParticle<ndim> &other = *neigh[i];
        const VectorND<ndim> eta = neigh[i]->disp - center->disp;
        const VectorND<ndim> xi = neigh[i]->coord - center->coord;
        const VectorND<ndim> zeta = xi + eta;
        Nmat += omega[i] * zeta.bun(xi) * center->vol[i];
    }

    MatrixND<ndim, ndim> Fmat = Nmat * Kinv;

    for (int i = 0; i < numfam; i++)
    {
        const VectorND<ndim> eta = neigh[i]->disp - center->disp;
        const VectorND<ndim> xi = neigh[i]->coord - center->coord;
        const VectorND<ndim> zeta = xi + eta;

        const MatrixND<ndim, ndim> Amat = this->get_A(xi);
        const MatrixND<ndim, ndim> Bmat = zeta.bun(xi) / xi.dot(xi);
        // Fmat <- Fmat * (I - xi\otimes xi / |xi|^2) + zeta\otimes xi / |xi|^2
        materials[i]->set_strain(Fmat * Amat + Bmat);
    }
}

template <int ndim, int maxfam>
MatrixND<ndim, ndim>
NosbProj<ndim, maxfam>::sum_PKinv()
{

    MatrixND<ndim, ndim> Qmat{{0.0}};
    for (int i = 0; i < numfam; i++)
    {
        const VectorND<ndim> xi = neigh[i]->coord - center->coord;
        const MatrixSD<ndim>& P = materials[i]->get_stress();
        const MatrixND<ndim, ndim> A = this->get_A(xi);
        // the correct formula for Qmat is 
        // Qmat += sum_1^numfam { w[i] * P * A * Kinv * vol[i] }
        // Since w[i] and omega[i] are the same in this formulation, we can use omega[i] instead
        Qmat.addMatrix(P * A * Kinv, omega[i] * center->vol[i]);


//      opserr  << Matrix(Qmat) << "\n";
    }
    return Qmat;
}

template <int ndim, int maxfam>
VectorND<ndim>
NosbProj<ndim, maxfam>::bond_force(int i, MatrixND<ndim, ndim> &Qmat)
{
    const VectorND<ndim> xi = neigh[i]->coord - center->coord;
    VectorND<ndim> T2 = this->get_T2(i, xi);

    // the correct formula for T[i] is
    // T[i] = omega[i] * Qmat * xi + w[i]*Pmat*xi/|xi|^2
    // Since w[i] and omega[i] are the same in this formulation and T2 = Pmat*xi/|xi|^2,
    // we can use the following formula instead:
    // T[i] = omega[i] * (Qmat * xi + T2)
    // Because the peridynamics force density at writes as
    // pforce = sum_1^numfam { T[i] * vol[i] }
    // T[i] is returned as 
    // T[i] = omega[i] * vol[i] * correction[i] * (Qmat * xi + T2)
    VectorND<ndim> T_j = Qmat * xi + T2;

    double w = omega[i] * center->vol[i] * center->correction[i];

    return T_j * w;
}
