#include <tcl.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <cstring>
#include <PeriDomain.h>
#include <PeriDomainBase.h>
// #include <typeinfo>

// Declare the function 'static' to convey that
// it is not expected to be used outside this file.
static int
Tcl_PeriInit(ClientData cd, Tcl_Interp* interp,
             int argc, const char** const argv)
{
  
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
  

  // clean up memory
  delete domain;

  return TCL_OK;
}


int
Tcl_Peri(ClientData cd, Tcl_Interp* interp,
         int argc, const char** const argv)
{
  // TODO Ensure argv[1] exists
  if (strcmp(argv[1], "init") == 0)
    return Tcl_PeriInit(cd, interp, argc, argv);
    
  return TCL_ERROR;
}
