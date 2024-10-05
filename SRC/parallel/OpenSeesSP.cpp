//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
//                                   OpenSeesSP                               
//
//===----------------------------------------------------------------------===//

extern "C" {
#include <tcl.h>
}


#include "commands.h"


#include <stdio.h>
#include <string.h>

#include <PartitionedDomain.h>
#include <MPI_MachineBroker.h>
#include <ShadowSubdomain.h>
#include <ActorSubdomain.h>
#include <TclPackageClassBroker.h>
#include <DomainPartitioner.h>

#include <mpi.h>

int Init_OpenSees(Tcl_Interp *interp);

void Init_MachineRuntime(Tcl_Interp* interp, MachineBroker* theMachineBroker);
void Init_PartitionRuntime(Tcl_Interp* interp, MachineBroker*, FEM_ObjectBroker*);

static int 
doNothing(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  return TCL_OK;
}


int
Openseessp_Init(Tcl_Interp* interp)
{
  int argc = 0; 
  char **argv = nullptr;

  MachineBroker* theMachineBroker = new MPI_MachineBroker(0, argc, argv);
  FEM_ObjectBroker* theBroker = new TclPackageClassBroker();
  theMachineBroker->setObjectBroker(theBroker);

  int pid = theMachineBroker->getPID();
  int np = theMachineBroker->getNP();

  //
  // depending on rank we do something
  //
  if (pid != 0) {

    //
    // on secondary processes we spin waiting to create & run actors
    //
    fprintf(stderr, "Secondary Process Running %d\n", pid);
    theMachineBroker->runActors();

  } else {

    //
    // on process 0 we create some ShadowSubdomains & then start the OpenSees interpreter
    //
    fprintf(stderr, "Primary Process Running OpenSees Interpreter %d\n", pid);   


    // Add machine commands (getPID, getNP, etc);
    Init_MachineRuntime(interp, theMachineBroker);
    Init_OpenSees(interp);
    Init_PartitionRuntime(interp, theMachineBroker, theBroker);

    Tcl_CreateCommand(interp, "barrier",   &doNothing, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

    // TODO
//  // some clean up to shut the remotes down if still running
//  theDomain.clearAll();
//  
//  // shutdown the remote machines
//  theMachineBroker->shutdown();
  }
  return 0;
}


