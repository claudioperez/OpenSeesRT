#include <tcl.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <cstring>
#include <PeriDomain.h>
#include <PeriDomainBase.h>
// #include <typeinfo>


Tcl_CmdProc Tcl_Peri;

// Declare the function 'static' to convey that
// it is not expected to be used outside this file.
static int
Tcl_PeriInit(ClientData cd, Tcl_Interp* interp,
             int argc, const char** const argv)
{
  // If init has already been called, clean up (free) the previous
  // domain object.
  if (cd != nullptr){
    PeriDomainBase *domain = static_cast<PeriDomainBase*>(cd);
    delete domain;
  }
  
  // Parse arguments and create a domain
  // Note that argv[0] holds the name of the command.
  
  // default values
  int ndim    = 0;
  int totnode = 0;
  int maxfam  = 0;
  char plane_type = 'x';
  PeriDomainBase *domain;

  int argi = 2;

  // if the size of argv is greater than 2, convert
  // the 3rd element of the array (ie argv[2]) to an integer
  // and store it in totnode
  if (argc > argi) {
    // If this returns TCL_ERROR it means that argv[2] couldnt be
    // parsed as an integer.
    if (Tcl_GetInt(interp, argv[argi], &ndim) == TCL_ERROR) {
      printf("ERROR in peri init: Couldnt parse argv[2] as an integer\n");
      return TCL_ERROR;
    }
    argi++;
  }

  // Now do the same for the 4th argument (argv[3])
  if (argc > argi) {
    if (Tcl_GetInt(interp, argv[argi], &totnode) == TCL_ERROR) {
      printf("ERROR in peri init: Couldnt parse argv[3] as an integer\n");
      return TCL_ERROR;
    }
    argi++;
  }

  // Now do the same for the 5th argument (argv[4])
  if (argc > argi) {
    if (Tcl_GetInt(interp, argv[argi], &maxfam) == TCL_ERROR) {
      printf("ERROR in peri init: Couldnt parse argv[4] as an integer\n");
      return TCL_ERROR;
    }
    argi++;
  }
  // argv[5] is only required for 2D
  if (argc > argi) {
    plane_type = char(argv[argi][0]);
  }

  printf("Creating %d D domain with totnode=%d and maxfam=%d...\n", ndim, totnode, maxfam);
  // Allocate a new domain
  if (ndim == 2) {
    if (plane_type != 's' && plane_type != 'e') {
      printf("ERROR in peri init: if ndim == 2, correct plane_type ('s' or 'e') is required.\n");
      return TCL_ERROR;
    }
    else {
      domain = new PeriDomain<2>(totnode, maxfam);
    }
    
  }
  else if (ndim == 3) {
    if (plane_type != 'x') {
      printf("EERROR in peri init: if ndim == 3, plane_type is not required.\n");
      return TCL_ERROR;
    }
    else {
      domain = new PeriDomain<3>(totnode, maxfam);
    }
  }
  else {
    printf("ERROR in peri init: ndim should be 2 or 3\n");
    return TCL_ERROR;
  }

  if (ndim == 3){
    printf("3D dmain with totnode=%d, maxfam=%d created.\n", domain->totnode, domain->maxfam);
  }else if (ndim == 2 && plane_type == 's'){
    printf("2D plane stress dmain with totnode=%d, maxfam=%d created.\n", domain->totnode, domain->maxfam);
  }else if (ndim == 2 && plane_type == 'e'){
    printf("2D plane strain dmain with totnode=%d, maxfam=%d created.\n", domain->totnode, domain->maxfam);
  }
  

  // Store the pointer to the domain in the ClientData
  Tcl_CreateCommand(interp, "peri", Tcl_Peri, static_cast<ClientData>(domain), nullptr);

  return TCL_OK;
}

static int
Tcl_PeriSetNode(ClientData cd, Tcl_Interp* interp,
                 int argc, const char** const argv)
{
  // decalre a pointer to a PeriDomainBase object
  PeriDomainBase* domain_base = static_cast<PeriDomainBase*>(cd);
  int     ind = 0, ndim = 0, argi = 2;

  ndim = domain_base->getNDM();
  if (ndim == 2) {
    // cast the pointer to a PeriDomain<2> object
    PeriDomain<2>* domain = static_cast<PeriDomain<2>*>(domain_base);
    std::array<double, 2> coord;

    // Parse arguments and set the coordinates of a particle
    // Note that argv[0] == "peri", argv[1] == "node"
    // Therefore, the first argument is argv[2]

    if (argc > argi) {
      // If this returns TCL_ERROR it means that argv[2] couldnt be
      // parsed as an integer.
      if (Tcl_GetInt(interp, argv[argi], &ind) == TCL_ERROR) {
        printf("ERROR in peri node: Couldnt parse argv[2] as the index\n");
        return TCL_ERROR;
      }
      argi++;
    }

    for (int j = 0; j < ndim; j++) {
      if (argc > argi) {
        if (Tcl_GetDouble(interp, argv[argi], &coord[j]) == TCL_ERROR) {
          printf("ERROR in peri node: Couldnt parse argv[%d] as a double\n", argi);
          return TCL_ERROR;
        }
        argi++;
      }else{
        printf("ERROR in peri node: Not enough arguments\n");
        return -1;
      }
    }
    // Set the coordinates of the particle at index i
    domain->set_coord(ind, coord);

  } else if (ndim == 3) {
    // cast the pointer to a PeriDomain<3> object
    PeriDomain<3>* domain = static_cast<PeriDomain<3>*>(domain_base);
    std::array<double, 3> coord;

    // Parse arguments and set the coordinates of a particle
    // Note that argv[0] == "peri", argv[1] == "node"
    // Therefore, the first argument is argv[2]

    if (argc > argi) {
      // If this returns TCL_ERROR it means that argv[2] couldnt be
      // parsed as an integer.
      if (Tcl_GetInt(interp, argv[argi], &ind) == TCL_ERROR) {
        printf("ERROR in peri node: Couldnt parse argv[2] as the index\n");
        return TCL_ERROR;
      }
      argi++;
    }

    for (int j = 0; j < ndim; j++) {
      if (argc > argi) {
        if (Tcl_GetDouble(interp, argv[argi], &coord[j]) == TCL_ERROR) {
          printf("ERROR in peri node: Couldnt parse argv[%d] as a double\n", argi);
          return TCL_ERROR;
        }
        argi++;
      }else{
        printf("ERROR in peri snode: Not enough arguments\n");
        return -1;  // QUESTION: Is this the correct return value?
      }
    }
    // Set the coordinates of the particle at index i
    domain->set_coord(ind, coord);
  }

  return TCL_OK;
}

static int
Tcl_PeriPrintNode(ClientData cd, Tcl_Interp* interp,
                 int argc, const char** const argv)
{
  // decalre a pointer to a PeriDomainBase object
  PeriDomainBase* domain_base = static_cast<PeriDomainBase*>(cd);
  int     ind = 0, ndim = 0, argi = 3;

  ndim = domain_base->getNDM();
  if (ndim == 2) {
    // cast the pointer to a PeriDomain<2> object
    PeriDomain<2>* domain = static_cast<PeriDomain<2>*>(domain_base);
    if (argc == 3) {
      // print the coordinates of all the particles
      for (int i = 0; i < domain->totnode; i++) {
        printf("Node %d: (%f, %f)\n", i, domain->pts[i].coord[0], domain->pts[i].coord[1]);
      }
    }else if (argc > argi) {
      // Parse arguments and set the coordinates of a particle
      // Note that argv[0] == "peri", argv[1] == "prin", argv[2] == "node"
      // Therefore, the first argument is argv[3]

      // If this returns TCL_ERROR it means that argv[3] couldnt be
      // parsed as an integer.
      if (Tcl_GetInt(interp, argv[argi], &ind) == TCL_ERROR) {
        printf("ERROR in peri prin node: Couldnt parse argv[3] as the index\n");
        return TCL_ERROR;
      }
      // Print the coordinates of the particle at index i
      printf("Node %d: (%f, %f)\n", ind, domain->pts[ind].coord[0], domain->pts[ind].coord[1]);
    }else{
      printf("ERROR in peri prin node: Not enough arguments\n");
      return -1;
    }
  } else if (ndim == 3) {
    // cast the pointer to a PeriDomain<3> object
    PeriDomain<3>* domain = static_cast<PeriDomain<3>*>(domain_base);

    if (argc == 3) {
      // print the coordinates of all the particles
      for (int i = 0; i < domain->totnode; i++) {
        printf("Node %d: (%f, %f, %f)\n", i, domain->pts[i].coord[0], domain->pts[i].coord[1], domain->pts[i].coord[2]);
      }
    }else if (argc > argi) {
      // Parse arguments and set the coordinates of a particle
      // Note that argv[0] == "peri", argv[1] == "prin", argv[2] == "node"
      // Therefore, the first argument is argv[3]

      // If this returns TCL_ERROR it means that argv[3] couldnt be
      // parsed as an integer.
      if (Tcl_GetInt(interp, argv[argi], &ind) == TCL_ERROR) {
        printf("ERROR in peri prin node: Couldnt parse argv[3] as the index\n");
        return TCL_ERROR;
      }
      // Print the coordinates of the particle at index i
      printf("Node %d: (%f, %f, %f)\n", ind, domain->pts[ind].coord[0], domain->pts[ind].coord[1], domain->pts[ind].coord[2]);
    }
  }
  return TCL_OK;
} 

int
Tcl_Peri(ClientData cd, Tcl_Interp* interp,
         int argc, const char** const argv)
{
  // TODO Ensure argv[1] exists
  if (strcmp(argv[1], "init") == 0){
    return Tcl_PeriInit(cd, interp, argc, argv);
  }
  
  if (cd == nullptr) {
    printf("ERROR: ClientData is null.\n");
    return -1;  // QUESTION: Is this the correct return value?
  }

  if (strcmp(argv[1], "node") == 0){
    return Tcl_PeriSetNode(cd, interp, argc, argv);
  }

  if (strcmp(argv[1], "prin") == 0){
    if (strcmp(argv[2], "node") == 0){
      return Tcl_PeriPrintNode(cd, interp, argc, argv);
    }
  }
    
    
  return TCL_OK;
}
