#include <tcl.h>
#include <stdio.h>

#include <ParticleDomain.h>

int
Tcl_PeridynamicsCommands(ClientData cd, Tcl_Interp* interp,
                         int argc, const char** const argv)
{

  // allocate a new domain
  ParticleDomain *domain = new ParticleDomain(1, 2);


  printf("Hello world\n");


  // Loop over the command argument array argv and print each one
  for (int i=0; i<argc; i++)
    printf("%s\n", argv[i]);


  domain->print(0);

  // clean up memory
  delete domain;

  return TCL_OK;
}
