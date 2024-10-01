#include <tcl.h>
#include <cmath>
#include <stdio.h>
#include <vector>
#include <string>
#include <cstring>
#include <Logging.h>
#include <PeriDomain.h>
#include <PeriElement.h>
#include <PeriDomainBase.h>

#include "mate/ElasticIsotropic.h"
#include <PeriParticle.h>
#include <NosbProj.h>

#include <threads/thread_pool.hpp>

Tcl_CmdProc Tcl_Peri;

// Declare the function 'static' to convey that
// it is not expected to be used outside this file.
static int
Tcl_PeriInit(ClientData cd, Tcl_Interp *interp,
             int argc, const char **const argv)
{
    // If init has already been called, clean up (free) the previous
    // domain object.
    if (cd != nullptr)
    {
        PeriDomainBase *domain = static_cast<PeriDomainBase *>(cd);
        delete domain;
    }

    // Parse arguments and create a domain
    // Note that argv[0] holds the name of the command.

    // default values
    int ndim = 0;
    int totnode = 0;
    int maxfam = 0;
    char plane_type = 'x';
    PeriDomainBase *domain;

    int argi = 2;

    // if the size of argv is greater than 2, convert
    // the 3rd element of the array (ie argv[2]) to an integer
    // and store it in totnode
    if (argc > argi)
    {
        // If this returns TCL_ERROR it means that argv[2] couldnt be
        // parsed as an integer.
        if (Tcl_GetInt(interp, argv[argi], &ndim) == TCL_ERROR)
        {
            printf("ERROR in peri init: Couldnt parse argv[2] as an integer\n");
            return TCL_ERROR;
        }
        argi++;
    }

    // Now do the same for the 4th argument (argv[3])
    if (argc > argi)
    {
        if (Tcl_GetInt(interp, argv[argi], &totnode) == TCL_ERROR)
        {
            printf("ERROR in peri init: Couldnt parse argv[3] as an integer\n");
            return TCL_ERROR;
        }
        argi++;
    }

    // Now do the same for the 5th argument (argv[4])
    if (argc > argi)
    {
        if (Tcl_GetInt(interp, argv[argi], &maxfam) == TCL_ERROR)
        {
            printf("ERROR in peri init: Couldnt parse argv[4] as an integer\n");
            return TCL_ERROR;
        }
        argi++;
    }
    // argv[5] is only required for 2D
    if (argc > argi)
    {
        plane_type = char(argv[argi][0]);
    }

    printf("Creating %d D domain with totnode=%d and maxfam=%d...\n", ndim, totnode, maxfam);
    // Allocate a new domain
    if (ndim == 2)
    {
        if (plane_type != 's' && plane_type != 'e')
        {
            printf("ERROR in peri init: if ndim == 2, correct plane_type ('s' or 'e') is required.\n");
            return TCL_ERROR;
        }
        else
        {
            domain = new PeriDomain<2>(totnode, maxfam);
        }
    }
    else if (ndim == 3)
    {
        if (plane_type != 'x')
        {
            printf("EERROR in peri init: if ndim == 3, plane_type is not required.\n");
            return TCL_ERROR;
        }
        else
        {
            domain = new PeriDomain<3>(totnode, maxfam);
        }
    }
    else
    {
        printf("ERROR in peri init: ndim should be 2 or 3\n");
        return TCL_ERROR;
    }

    if (ndim == 3)
    {
        printf("3D dmain with totnode=%d, maxfam=%d created.\n", domain->totnode, domain->maxfam);
    }
    else if (ndim == 2 && plane_type == 's')
    {
        printf("2D plane stress dmain with totnode=%d, maxfam=%d created.\n", domain->totnode, domain->maxfam);
    }
    else if (ndim == 2 && plane_type == 'e')
    {
        printf("2D plane strain dmain with totnode=%d, maxfam=%d created.\n", domain->totnode, domain->maxfam);
    }

    // Store the pointer to the domain in the ClientData
    Tcl_CreateCommand(interp, "peri", Tcl_Peri, static_cast<ClientData>(domain), nullptr);

    return TCL_OK;
}



template <int ndim>
static int
Tcl_PeriSetNode(PeriDomain<ndim> *domain, Tcl_Interp *interp,
                int argc, const char **const argv)
{
    int ind = 0, argi = 2;
    std::array<double, ndim> coord;

    // Parse arguments and set the coordinates of a particle
    // Note that argv[0] == "peri", argv[1] == "node"
    // Therefore, the first argument is argv[2]

    if (argc > argi)
    {
        // If this returns TCL_ERROR it means that argv[2] couldnt be
        // parsed as an integer.
        if (Tcl_GetInt(interp, argv[argi], &ind) == TCL_ERROR)
        {
            printf("ERROR in peri node: Couldnt parse argv[2] as the index\n");
            return TCL_ERROR;
        }
        argi++;
    }

    for (int j = 0; j < ndim; j++)
    {
        if (argc > argi)
        {
            if (Tcl_GetDouble(interp, argv[argi], &coord[j]) == TCL_ERROR)
            {
                printf("ERROR in peri node: Couldnt parse argv[%d] as a double\n", argi);
                return TCL_ERROR;
            }
            argi++;
        }
        else
        {
            printf("ERROR in peri node: Not enough arguments\n");
            return -1; // **QUESTION: Is this the correct return value?**
        }
    }
    // Set the coordinates of the particle at index i
    domain->set_coord(ind, coord);

    return TCL_OK;
}

template <int ndim>
static int
Tcl_PeriPrintNode(PeriDomain<ndim> *domain, Tcl_Interp *interp,
                  int argc, const char **const argv)
{
    int ind = 0, argi = 3;

    if (argc == 3)
    {
        // print the coordinates of all the particles
        for (int i = 0; i < domain->totnode; i++)
        {
            printf("Node %d: (", i);
            for (int j = 0; j < ndim; j++)
            {
                printf("%f", domain->pts[i].coord[j]);
                if (j < ndim - 1)
                {
                    printf(", ");
                }
            }
            printf(")\n");
        }
    }
    else if (argc > argi)
    {
        // Parse arguments and set the coordinates of a particle
        // Note that argv[0] == "peri", argv[1] == "prin", argv[2] == "node"
        // Therefore, the first argument is argv[3]
        for (int i = argi; i < argc; i++)
        {
            // If this returns TCL_ERROR it means that argv[3] couldnt be
            // parsed as an integer.
            if (Tcl_GetInt(interp, argv[i], &ind) == TCL_ERROR)
            {
                printf("ERROR in peri prin node: Couldnt parse argv[3] as the index\n");
                return TCL_ERROR;
            }
            // Print the coordinates of the particle at index i
            printf("Node %d: (", ind);
            for (int j = 0; j < ndim; j++)
            {
                printf("%f", domain->pts[ind].coord[j]);
                if (j < ndim - 1)
                {
                    printf(", ");
                }
            }
            printf(")\n");
        }
    }
    else
    {
        printf("ERROR in peri prin node: Not enough arguments\n");
        return -1; // **QUESTION: Is this the correct return value?**
    }
    return TCL_OK;
}

template <int ndim>
static int
Tcl_PeriCreateFam(PeriDomain<ndim> *domain, Tcl_Interp *interp,
                  int argc, const char **const argv)
{
    int argi = 2;
    double delta = 0.0;

    // Parse arguments and set the coordinates of a particle
    // Note that argv[0] == "peri", argv[1] == "fam"
    // Therefore, the first argument is argv[2]

    if (argc > argi)
    {
        // If this returns TCL_ERROR it means that argv[2] couldnt be
        // parsed as a double.
        if (Tcl_GetDouble(interp, argv[argi], &delta) == TCL_ERROR)
        {
            printf("ERROR in peri fam: Couldnt parse argv[2] as a double\n");
            return TCL_ERROR;
        }
        argi++;
    }

    // Create families for each particle
    domain->create_fam(delta);

    return TCL_OK;
}


template <int ndim>
static int
Tcl_PeriSetVols(PeriDomain<ndim> *domain, Tcl_Interp *interp,
                int argc, const char **const argv)
{
    int ind = 0, argi = 2;
    double vol = 0.0;

    // Parse arguments and set the volume of a particle
    // Note that argv[0] == "peri", argv[1] == "vol"
    // Therefore, the first argument is argv[2]

    if (argc > argi)
    {
        // If this returns TCL_ERROR it means that argv[2] couldnt be
        // parsed as an integer.
        if (Tcl_GetInt(interp, argv[argi], &ind) == TCL_ERROR)
        {
            printf("ERROR in peri vol: Couldnt parse argv[2] as the index\n");
            return TCL_ERROR;
        }
        argi++;
    }

    if (argc > argi)
    {
        // If this returns TCL_ERROR it means that argv[3] couldnt be
        // parsed as a double.
        if (Tcl_GetDouble(interp, argv[argi], &vol) == TCL_ERROR)
        {
            printf("ERROR in peri vol: Couldnt parse argv[3] as a double\n");
            return TCL_ERROR;
        }
        argi++;
    }

    // Set the volume of the particle at index i
    domain->set_vols(ind, vol);

    return TCL_OK;
}

template <int ndim>
static int
Tcl_PeriPrintFam(PeriDomain<ndim> *domain, Tcl_Interp *interp,
                 int argc, const char **const argv)
{
    int ind = 0, argi = 3;

    if (argc == 3)
    {
        // print the number of families and their indices of all the particles
        for (int i = 0; i < domain->totnode; i++)
        {
            printf("Node %d: %d neighbors\n", i, domain->pts[i].numfam);
            for (int j = 0; j < domain->pts[i].numfam; j++)
            {
                printf("%d ", domain->pts[i].nodefam[j]);
            }
            printf("\n");
        }
    }
    else if (argc > argi)
    {
        // Parse arguments and set the coordinates of a particle
        // Note that argv[0] == "peri", argv[1] == "prin", argv[2] == "fam"
        // Therefore, the first argument is argv[3]
        for (int i = argi; i < argc; i++)
        {
            // If this returns TCL_ERROR it means that argv[3] couldnt be
            // parsed as an integer.
            if (Tcl_GetInt(interp, argv[i], &ind) == TCL_ERROR)
            {
                printf("ERROR in peri prin fam: Couldnt parse argv[3] as the index\n");
                return TCL_ERROR;
            }
            // Print the families of the particle at index i
            printf("Node %d: %d neighbors\n", ind, domain->pts[ind].numfam);
            for (int j = 0; j < domain->pts[ind].numfam; j++)
            {
                printf("%d ", domain->pts[ind].nodefam[j]);
            }
            printf("\n");
        }
    }
    else
    {
        printf("ERROR in peri prin fam: Not enough arguments\n");
        return -1; // **QUESTION: Is this the correct return value?**
    }
    return TCL_OK;
}


template <int ndim>
static int
Tcl_PeriCalcVols(PeriDomain<ndim> *domain, Tcl_Interp *interp,
                 int argc, const char **const argv)
{
    int argi = 2;
    double space = 0.0;

    // Parse arguments and set the volume of a particle
    // Note that argv[0] == "peri", argv[1] == "cvol"
    // Therefore, the first argument is argv[2]

    if (argc > argi)
    {
        // If this returns TCL_ERROR it means that argv[2] couldnt be
        // parsed as a double.
        if (Tcl_GetDouble(interp, argv[argi], &space) == TCL_ERROR)
        {
            printf("ERROR in peri cvol: Couldnt parse argv[2] as a double\n");
            return TCL_ERROR;
        }
        argi++;
    }

    // Calculate the volume of the horizons
    domain->calc_vols(space);

    return TCL_OK;
}


template<int ndim>
static int
Tcl_PeriElem(ClientData clientData, Tcl_Interp *interp,
                   int argc, const char **const argv)
{
//PeriDomain<ndim> *pdomain;
    PeriDomainBase *domain = static_cast<PeriDomainBase *>(clientData);
  PeriElement(1, pdomain);

  return TCL_OK;
}


template <int ndim>
static int
Tcl_PeriPrintVol(PeriDomain<ndim> *domain, Tcl_Interp *interp,
                 int argc, const char **const argv)
{
    int argi = 3, idx = 0;

    if (argc == 3)
    {
        // print the volume of all the particles
        for (int i = 0; i < domain->totnode; i++)
        {
            printf("Node %d: %f\n", i, domain->pts[i].vol_h);
            for (int ind = 0; ind < domain->pts[i].numfam; ind++)
            {
                printf("%7.2E ", domain->pts[i].vol[ind]);
            }
            printf("\n");
        }
    }
    else if (argc > argi)
    {
        // Parse arguments and set the coordinates of a particle
        // Note that argv[0] == "peri", argv[1] == "prin", argv[2] == "vol"
        // Therefore, the first argument is argv[3]
        for (int k = argi; k < argc; k++)
        {
            // If this returns TCL_ERROR it means that argv[3] couldnt be
            // parsed as an integer.
            if (Tcl_GetInt(interp, argv[k], &idx) == TCL_ERROR)
            {
                printf("ERROR in peri prin vol: Couldnt parse argv[3] as the index\n");
                return TCL_ERROR;
            }
            // Print the volume of the particle at index idx
            printf("Node %d: %f\n", idx, domain->pts[idx].vol_h);
            for (int ind = 0; ind < domain->pts[idx].numfam; ind++)
            {
                printf("%7.2E ", domain->pts[idx].vol[ind]);
            }
            printf("\n");
        }
    }
    else
    {
        printf("ERROR in peri prin vol: Not enough arguments\n");
        return -1; // **QUESTION: Is this the correct return value?**
    }
    return TCL_OK;
}


int
Tcl_PeriFormThreads(PeriDomain<3>& domain, Tcl_Interp* interp, int argc, const char** const argv)
{
  OpenSees::thread_pool threads{6};
  std::mutex resp_mutex;

  constexpr int ndim = 3;
  constexpr int maxfam = 1024;

  std::vector<NosbProj<3,maxfam>> nodefam;

  // Create families for specific NOSB type
  for (PeriParticle<ndim>& particle : domain.pts) {
    nodefam.emplace_back(&particle, domain, new ElasticIsotropic<ndim>(29e3, 0.2));
  }


  // Initialize shape tensor
  threads.submit_loop<unsigned int>(0, nodefam.size(), [&](int i){
    nodefam[i].init_shape();
  }).wait();


  // Form deformation gradients for trial
  threads.submit_loop<unsigned int>(0, nodefam.size(), [&](int i){
    nodefam[i].form_trial();
  }).wait();

  // Form force
  threads.submit_loop<unsigned int>(0, nodefam.size(), [&](int i){
    MatrixND<ndim,ndim> Q = nodefam[i].sum_PKinv();

    for (int j=0; j < nodefam[i].numfam; j++) {
      const VectorND<ndim> T_j = nodefam[i].bond_force(j, Q);

      const std::lock_guard<std::mutex> lock(resp_mutex);
      nodefam[i].center->pforce   += T_j;
      nodefam[i].neigh[j]->pforce -= T_j;
    }

  }).wait();

  //
  // TODO: Update disp
  //

  // Create and return force vector
  Tcl_Obj* list = Tcl_NewListObj(nodefam.size()*ndim, nullptr);
  for (PeriParticle<ndim>& particle : domain.pts) 
    for (int j = 0; j<ndim; j++)
      Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj(particle.pforce[j]));

  Tcl_SetObjResult(interp, list);
  return TCL_OK;
}


int
Tcl_PeriStep(PeriDomain<3>& domain, Tcl_Interp* interp, int argc, const char** const argv)
{
    constexpr int ndim = 3;
    const int nn = domain.pts.size();

    double eps_val = 1e-8;


    const int tt = 0;
    const double dt_hs = 0.1,
          damp_ratio = 0.1;

    Matrix velhalf(nn,ndim),
           velhalfold(nn, ndim),
           pforceold(nn, ndim);
    Vector dens(nn);



    double cdamp = 0.0, cdamp1 = 0.0, cdamp2 = 0.0;

    // calculate damping coefficient cdamp
    for (int i=0; i<nn; i++) {
        PeriParticle<ndim>& particle = domain.pts[i];
        for(int j = 0; j < ndim; j++) {
            if (abs(velhalfold(i, j)) > eps_val)  {
                double Kii = (particle.pforce[j] - pforceold(i, j)) / (velhalfold(i, j)*dt_hs);
                cdamp1 += particle.disp[j] * 1.0/dens(i) * (-Kii) * particle.disp[j];
            }
            cdamp2 = cdamp2 + particle.disp[j] * particle.disp[j];
        }
    }
    if (cdamp2 > eps_val)  {
        if ((cdamp1 / cdamp2) > 0.0)
            cdamp = 2.0 * sqrt(cdamp1 / cdamp2);
        else
            cdamp = 0.0;

    } else {
        cdamp = 0.0;
    }

    cdamp = std::min(cdamp, 1.9) * damp_ratio;

    for (int i=0; i<nn; i++) {
        PeriParticle<ndim>& particle = domain.pts[i];
        for(int j = 0; j < ndim; j++) {
            if (particle.is_disp_bound[j] == 0)  {
                if (tt == 1)  {
                      velhalf(i,j) = particle.vel[j] + 0.5 * dt_hs * (particle.pforce[j] + particle.bforce[j]) / dens(i);

                } else {
                      velhalf(i,j) = ((2.0 - cdamp * dt_hs) * velhalfold(i,j) +
                                       2.0 * dt_hs * (particle.pforce[j] + particle.bforce[j]) / dens(i))
                                     / (2.0 + cdamp * dt_hs);
                }
                particle.vel[j] = 0.5 * (velhalfold(i,j) + velhalf(i,j));
                particle.disp[j] = particle.disp[j] + velhalf(i,j) * dt_hs;
                velhalfold(i,j) = velhalf(i,j);
            } else  {
                particle.vel[j]   = particle.bdisp[j];
                velhalf(i, j)     = particle.bdisp[j];
                velhalfold(i, j)  = particle.bdisp[j];
                particle.disp[j] += particle.bdisp[j] * dt_hs;
            }
            pforceold(i, j) = particle.pforce[j];
        }
    }

}

int
Tcl_PeriForm(PeriDomain<3>& domain, Tcl_Interp* interp, int argc, const char** const argv)
{
  constexpr int ndim = 3;
  constexpr int maxfam = 1024;

  std::vector<NosbProj<ndim,maxfam>> nodefam;

  // Create families for specific NOSB type
  for (PeriParticle<3>& particle : domain.pts) {
    nodefam.emplace_back(&particle, domain, new ElasticIsotropic<ndim>(29e3, 0.2));
  }

  // Initialize shape tensor
  for (NosbProj<ndim,maxfam>& fam_i : nodefam)
    fam_i.init_shape();

  // Form deformation gradients for trial
  for (NosbProj<ndim,maxfam>& fam_i : nodefam)
    fam_i.form_trial();


  // Form force
  for (NosbProj<3,maxfam>& fam_i : nodefam) {

    MatrixND<ndim,ndim> Q = fam_i.sum_PKinv();

    for (int j=0; j < fam_i.numfam; j++) {
      const VectorND<ndim> T_j = fam_i.bond_force(j, Q);

      fam_i.center->pforce   += T_j;
      fam_i.neigh[j]->pforce -= T_j;

    }
  }

  //
  // TODO: Update disp
  //

  // Create and return force vector
  Tcl_Obj* list = Tcl_NewListObj(nodefam.size()*ndim, nullptr);
  for (PeriParticle<ndim>& particle : domain.pts) 
    for (int j = 0; j<ndim; j++)
      Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj(particle.pforce[j]));

  Tcl_SetObjResult(interp, list);
  return TCL_OK;
}

int Tcl_Peri(ClientData cd, Tcl_Interp *interp,
             int argc, const char **const argv)
{
    // TODO Ensure argv[1] exists
    if (strcmp(argv[1], "init") == 0)
    {
        return Tcl_PeriInit(cd, interp, argc, argv);
    }

    if (cd == nullptr)
    {
        printf("ERROR: ClientData is null.\n");
        return -1; // **QUESTION: Is this the correct return value?**
    }

    PeriDomainBase *domain_base = static_cast<PeriDomainBase *>(cd);

    if (strcmp(argv[1], "node") == 0)
    {
        int ndim = domain_base->getNDM();
        if (ndim == 2)
        {
            // cast the pointer to a PeriDomain<2> object
            PeriDomain<2> *domain = static_cast<PeriDomain<2> *>(domain_base);
            return Tcl_PeriSetNode(domain, interp, argc, argv);
        }
        else if (ndim == 3)
        {
            // cast the pointer to a PeriDomain<3> object
            PeriDomain<3> *domain = static_cast<PeriDomain<3> *>(domain_base);
            return Tcl_PeriSetNode(domain, interp, argc, argv);
        }
    }

    if (strcmp(argv[1], "fam") == 0)
    {
        int ndim = domain_base->getNDM();
        if (ndim == 2)
        {
            // cast the pointer to a PeriDomain<2> object
            PeriDomain<2> *domain = static_cast<PeriDomain<2> *>(domain_base);
            return Tcl_PeriCreateFam(domain, interp, argc, argv);
        }
        else if (ndim == 3)
        {
            // cast the pointer to a PeriDomain<3> object
            PeriDomain<3> *domain = static_cast<PeriDomain<3> *>(domain_base);
            return Tcl_PeriCreateFam(domain, interp, argc, argv);
        }
    }

    if (strcmp(argv[1], "svol") == 0)
    {
        int ndim = domain_base->getNDM();
        if (ndim == 2)
        {
            // cast the pointer to a PeriDomain<2> object
            PeriDomain<2> *domain = static_cast<PeriDomain<2> *>(domain_base);
            return Tcl_PeriSetVols(domain, interp, argc, argv);
        }
        else if (ndim == 3)
        {
            // cast the pointer to a PeriDomain<3> object
            PeriDomain<3> *domain = static_cast<PeriDomain<3> *>(domain_base);
            return Tcl_PeriSetVols(domain, interp, argc, argv);
        }
    }

    if (strcmp(argv[1], "cvol") == 0)
    {
        int ndim = domain_base->getNDM();
        if (ndim == 2)
        {
            // cast the pointer to a PeriDomain<2> object
            PeriDomain<2> *domain = static_cast<PeriDomain<2> *>(domain_base);
            return Tcl_PeriCalcVols(domain, interp, argc, argv);
        }
        else if (ndim == 3)
        {
            // cast the pointer to a PeriDomain<3> object
            PeriDomain<3> *domain = static_cast<PeriDomain<3> *>(domain_base);
            return Tcl_PeriCalcVols(domain, interp, argc, argv);
        }
    }

    if (strcmp(argv[1], "form") == 0)
    {
        int ndim = domain_base->getNDM();
        if (ndim == 3)
        {
            // cast the pointer to a PeriDomain<3> object
            PeriDomain<3> *domain = static_cast<PeriDomain<3> *>(domain_base);
      if (argc > 2)
        return Tcl_PeriFormThreads(*domain, interp, argc, argv);
      else
        return Tcl_PeriForm(*domain, interp, argc, argv);
    }

  }

    if (strcmp(argv[1], "prin") == 0)
    {

        int ndim = domain_base->getNDM();
        if (ndim == 2)
        {
            // cast the pointer to a PeriDomain<2> object
            PeriDomain<2> *domain = static_cast<PeriDomain<2> *>(domain_base);
            if (strcmp(argv[2], "node") == 0)
            {
                return Tcl_PeriPrintNode(domain, interp, argc, argv);
            }
            else if (strcmp(argv[2], "fam") == 0)
            {
                return Tcl_PeriPrintFam(domain, interp, argc, argv);
            }
            else if (strcmp(argv[2], "vol") == 0)
            {
                return Tcl_PeriPrintVol(domain, interp, argc, argv);
            }
        }
        else if (ndim == 3)
        {
            // cast the pointer to a PeriDomain<3> object
            PeriDomain<3> *domain = static_cast<PeriDomain<3> *>(domain_base);
            if (strcmp(argv[2], "node") == 0)
            {
                return Tcl_PeriPrintNode(domain, interp, argc, argv);
            }
            else if (strcmp(argv[2], "fam") == 0)
            {
                return Tcl_PeriPrintFam(domain, interp, argc, argv);
            }
            else if (strcmp(argv[2], "vol") == 0)
            {
                return Tcl_PeriPrintVol(domain, interp, argc, argv);
            }
        }
    }
    return TCL_OK;
}
