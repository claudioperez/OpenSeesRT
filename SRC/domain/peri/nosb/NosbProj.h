#include <Mate.h>
#include <NosbBase.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <Vector3D.h>
#include <Matrix3D.h>

using namespace OpenSees;

template <int ndim, int maxfam>
class NosbProj : public NosbBase<ndim>
{
public:
    int numfam;
    double scale, delta;
    MatrixND<ndim, ndim> Kinv; // note that Kmat is invariant in this formulation, but may not be in general
    PeriParticle<ndim> *center;
    std::array<PeriParticle<ndim> *, maxfam> neigh;
    std::array<Mate<ndim> *, maxfam> materials;
    std::array<double, maxfam> omega;

public:
    NosbProj(PeriParticle<ndim> *center, PeriDomain<ndim> &domain, Mate<ndim> *material);

    virtual ~NosbProj() {};

    virtual void init_shape() override final;

    virtual void form_trial() final;

    virtual MatrixND<ndim, ndim> sum_PKinv() final;

    virtual VectorND<ndim> bond_force(int i, MatrixND<ndim, ndim> &Qmat) final;

    // double bond_omega(int i);
    void bond_omega(int i, int w_type);

private:
    MatrixND<ndim, ndim> get_A(const VectorND<ndim> &xi);

    MatrixND<ndim, ndim> get_B(const VectorND<ndim> &xi, const VectorND<ndim> &zeta);

    VectorND<ndim> get_T2(int i, const VectorND<ndim> &xi);
};

#include <NosbProj.tpp> // Include the implementation of the template class
