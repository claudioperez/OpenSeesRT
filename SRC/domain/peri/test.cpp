#include <PeriParticle.h>
#include <PeriDomain.h>


#include <NosbProj.h>

int
test(PeriDomain<3>& domain)
{
  constexpr int ndim = 3;

  std::vector<NosbProj<3,10>> nodefam;

  // Create families for specific NOSB type
  for (PeriParticle<3>& particle : domain.pts) {
    nodefam.emplace_back(&particle, domain, new Mate());
  }

  // Initialize shape tensor
  for (NosbBase<3>& fam_i : nodefam)
    fam_i.init_shape();

  // Form deformation gradients for trial
  for (NosbBase<3>& fam_i : nodefam)
    fam_i.form_trial();


  for (NosbProj<3,10>& fam_i : nodefam) {

    MatrixND<ndim,ndim> Q = fam_i.sum_PKinv();

    for (int j=0; j < fam_i.numfam; j++) {
      const VectorND<ndim> T_j = fam_i.bond_force(j, Q);

      fam_i.center->pforce     += T_j;
      fam_i.neigh[j]->pforce -= T_j;

    }
  }

  //
  return 0;
}

