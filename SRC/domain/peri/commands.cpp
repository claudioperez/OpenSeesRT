#include <tcl.h>
#include <stdio.h>
#include <PeriDomain.h>
// #include <typeinfo>

int
Tcl_Peri(ClientData cd, Tcl_Interp* interp,
         int argc, const char** const argv)
{

  // //
  // printf("Hello world\n");

  // // Loop over the command argument array argv and print each one
  // for (int i=0; i<argc; i++)
  //   printf("%s\n", argv[i]);




  //
  // Parse arguments and create a domain
  //
  // Note that argv[0] holds the name of the command.
  
  // default values
  int ndim    = 0;
  int totnode = 0;
  int maxfam  = 0;
  char plane_type = 'x';
  PeriDomain *domain;

  // if the size of argv is greater than 1, convert
  // the second element of the array (ie argv[1]) to an integer
  // and store it in totnode
  if (argc > 1) {
    // If this returns TCL_ERROR it means that argv[1] couldnt be
    // parsed as an integer.
    if (Tcl_GetInt(interp, argv[1], &ndim) == TCL_ERROR) {
      printf("ERROR in peri: Couldnt parse argv[1] as an integer\n");
      return -1;
    }
  }

  // Now do the same for the third argument (argv[2])
  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &totnode) == TCL_ERROR) {
      printf("ERROR in peri: Couldnt parse argv[2] as an integer\n");
      return -1;
    }
  }

  // Now do the same for the fourth argument (argv[3])
  if (argc > 3) {
    if (Tcl_GetInt(interp, argv[3], &maxfam) == TCL_ERROR) {
      printf("ERROR in peri: Couldnt parse argv[3] as an integer\n");
      return -1;
    }
  }
  // argv[4] is optional
  if (argc > 4) {
    plane_type = char(argv[4][0]);
    // printf("Type is %s, value is %c\n", typeid(plane_type).name(), plane_type);
  }

  printf("Creating %d D domain with totnode=%d and maxfam=%d...\n", ndim, totnode, maxfam);
  // Allocate a new domain
  if (ndim == 2) {
    if (plane_type != 's' && plane_type != 'e') {
      printf("ERROR in peri: if ndim == 2, correct plane_type ('s' or 'e') is required.\n");
      return -1;
    }
    else {
      domain = new PeriDomain(ndim, totnode, maxfam, plane_type);
    }
    
  }
  else if (ndim == 3) {
    if (plane_type != 'x') {
      printf("EERROR in peri: if ndim == 3, plane_type is not required.\n");
      return -1;
    }
    else {
      domain = new PeriDomain(ndim, totnode, maxfam);
    }
  }
  else {
    printf("ERROR in peri: ndim should be 2 or 3\n");
    return -1;
  }

  // Print the domain
  // domain->print(0);
  printf("%d D dmain with totnode=%d, maxfam=%d created.\n", domain->ndim, domain->totnode, domain->maxfam);

  // clean up memory
  delete domain;

  return TCL_OK;
}
