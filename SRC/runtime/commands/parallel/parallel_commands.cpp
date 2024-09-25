//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
//
#include <tcl.h>
#include <OPS_Globals.h>
#include <mpi.h>
#include <MachineBroker.h>
#include <MPI_MachineBroker.h>

// #define _PARALLEL_SP
#define _PARALLEL_MP

#ifdef _PARALLEL_SP
// #  include <DistributedDisplacementControl.h>
// #  include <ShedHeaviest.h>
// #  include <MPIDiagonalSOE.h>
// #  include <MPIDiagonalSolver.h>
#  include <ShadowSubdomain.h>
#  include <Metis.h>
#  include <DomainPartitioner.h>
#  include <domain/domain/partitioned/PartitionedDomain.h>
#  include <GraphPartitioner.h>
#  include <TclPackageClassBroker.h>
#  include <Subdomain.h>
#  include <SubdomainIter.h>
#  include <MachineBroker.h>
#  include <StaticDomainDecompositionAnalysis.h>
#  include <TransientDomainDecompositionAnalysis.h>

// #  define MPIPP_H
// #  include <DistributedSuperLU.h>
// #  include <DistributedProfileSPDLinSOE.h>

   struct PartitionRuntime {
     int  PARALLEL_SP = 0;
     int  num_subdomains = 0;
     bool PARTITIONED = false;
     bool USING_MAIN_DOMAIN = false;
     bool setMPIDSOEFlag = false;
     int  MAIN_DOMAIN_PARTITION_ID = 0;
     PartitionedDomain     theDomain;
     DomainPartitioner     *DOMAIN_PARTITIONER = nullptr;
     GraphPartitioner      *GRAPH_PARTITIONER  = nullptr;
     LoadBalancer          *BALANCER           = nullptr;
     TclPackageClassBroker *broker             = nullptr;
     MachineBroker         *machine            = nullptr;
     Channel               **theChannels       = nullptr;  
   };

// MachineBroker *theMachineBroker = 0;
   int  OPS_PARALLEL_SP = 0;
   bool OPS_PARTITIONED = false;
   bool OPS_USING_MAIN_DOMAIN = false;
   bool setMPIDSOEFlag = false;
   int  OPS_MAIN_DOMAIN_PARTITION_ID = 0;
   PartitionedDomain     theDomain;
   DomainPartitioner     *OPS_DOMAIN_PARTITIONER = nullptr;
   GraphPartitioner      *OPS_GRAPH_PARTITIONER  = nullptr;
   LoadBalancer          *OPS_BALANCER           = nullptr;
   TclPackageClassBroker *OPS_OBJECT_BROKER      = nullptr;
   MachineBroker         *OPS_MACHINE            = nullptr;
   Channel               **OPS_theChannels       = nullptr;  

#  elif defined(_PARALLEL_MP)

  bool setMPIDSOEFlag = false;
  
  // parallel analysis
  
// Domain theDomain;
#endif

int getPID(ClientData,  Tcl_Interp *, int, TCL_Char ** const argv);
int getNP( ClientData,  Tcl_Interp *, int, TCL_Char ** const argv);
int opsBarrier(ClientData, Tcl_Interp *, int, TCL_Char ** const argv);
int opsSend(ClientData, Tcl_Interp *, int, TCL_Char ** const argv);
int opsRecv(ClientData, Tcl_Interp *, int,TCL_Char ** const argv);
int opsPartition(ClientData, Tcl_Interp *, int, TCL_Char ** const argv);
int wipePP(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

void Init_Parallel(Tcl_Interp* interp)
{
  MachineBroker* theMachineBroker = new MPI_MachineBroker(nullptr, 0, nullptr);
  Tcl_CreateCommand(interp, "getNP",     &getNP,   (ClientData)theMachineBroker, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "getPID",    &getPID,  (ClientData)theMachineBroker, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "send",      &opsSend, (ClientData)theMachineBroker, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "recv",      &opsRecv, (ClientData)theMachineBroker, (Tcl_CmdDeleteProc *)NULL);

#ifdef _PARALLEL_SP
  PartitionRuntime *part = new PartitionRuntime();
  Tcl_CreateCommand(interp, "barrier",   &opsBarrier,   (ClientData)part, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "partition", &opsPartition, (ClientData)part, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "wipePP",    &wipePP,       (ClientData)part, (Tcl_CmdDeleteProc *)NULL);
#else
  Tcl_CreateCommand(interp, "barrier",   &opsBarrier, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "partition", &opsPartition, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
#endif
}


int
getPID(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  int pid = 0;

  MachineBroker* theMachineBroker = (MachineBroker*)clientData;

  if (theMachineBroker != nullptr)
    pid = theMachineBroker->getPID();

  // now we copy the value to the tcl string that is returned
  char buffer[30];
  sprintf(buffer, "%d", pid);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
getNP(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  int np = 1;
  MachineBroker* theMachineBroker = (MachineBroker*)clientData;

  if (theMachineBroker != nullptr)
    np = theMachineBroker->getNP();

  // now we copy the value to the tcl string that is returned
  char buffer[30];
  sprintf(buffer, "%d", np);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}


#ifdef _PARALLEL_SP
static int
partitionModel(PartitionRuntime& part, int eleTag)
{
  if (part.PARTITIONED == true)
    return 0;

  int result = 0;

  if (part.theChannels != nullptr)
    delete[] part.theChannels;

  part.theChannels = new Channel *[part.num_subdomains];

  // create some subdomains
  for (int i = 1; i <= part.num_subdomains; i++) {
    if (i != OPS_MAIN_DOMAIN_PARTITION_ID) {
      ShadowSubdomain *theSubdomain =
          new ShadowSubdomain(i, *part.machine, *part.broker);
      part.theDomain.addSubdomain(theSubdomain);
      OPS_theChannels[i - 1] = theSubdomain->getChannelPtr();
    }
  }

  // create a partitioner & partition the domain
  if (part.DOMAIN_PARTITIONER == nullptr) {
    //      OPS_BALANCER = new ShedHeaviest();
    part.GRAPH_PARTITIONER = new Metis;
    // OPS_DOMAIN_PARTITIONER = new DomainPartitioner(*OPS_GRAPH_PARTITIONER, *OPS_BALANCER);
    part.DOMAIN_PARTITIONER = new DomainPartitioner(*part.GRAPH_PARTITIONER);
    part.theDomain.setPartitioner(part.DOMAIN_PARTITIONER);
  }

  result = part.theDomain.partition(part.num_subdomains, OPS_USING_MAIN_DOMAIN,
                                    OPS_MAIN_DOMAIN_PARTITION_ID, eleTag);

  if (result < 0)
    return result;

  OPS_PARTITIONED = true;

  DomainDecompositionAnalysis *theSubAnalysis;
  SubdomainIter &theSubdomains = part.theDomain.getSubdomains();
  Subdomain *theSub = nullptr;

  void* the_static_analysis = nullptr;

  // create the appropriate domain decomposition analysis
  while ((theSub = theSubdomains()) != nullptr) {
    if (the_static_analysis != nullptr) {
      theSubAnalysis = new StaticDomainDecompositionAnalysis(
          *theSub, *theHandler, *theNumberer, *the_analysis_model, *theAlgorithm,
          *theSOE, *theStaticIntegrator, theTest, false);

    } else {
      theSubAnalysis = new TransientDomainDecompositionAnalysis(
          *theSub, *theHandler, *theNumberer, *the_analysis_model, *theAlgorithm,
          *theSOE, *theTransientIntegrator, theTest, false);
    }
    theSub->setDomainDecompAnalysis(*theSubAnalysis);
  }

  return result;
}
#endif


int
opsBarrier(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
#ifdef _PARALLEL_MP
  return MPI_Barrier(MPI_COMM_WORLD);
#endif
  return TCL_OK;
}

int
opsSend(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  if (argc < 2)
    return TCL_OK;

  int otherPID = -1;
  MachineBroker* theMachineBroker = (MachineBroker*)clientData;
  int myPID = theMachineBroker->getPID();
  int np = theMachineBroker->getNP();
  const char *dataToSend = argv[argc - 1];
  int msgLength = strlen(dataToSend) + 1;

  const char *gMsg = dataToSend;
  //  strcpy(gMsg, dataToSend);

  if (strcmp(argv[1], "-pid") == 0 && argc > 3) {

    if (Tcl_GetInt(interp, argv[2], &otherPID) != TCL_OK) {
      opserr << "send -pid pid? data? - pid: " << argv[2] << " invalid\n";
      return TCL_ERROR;
    }

    if (otherPID > -1 && otherPID != myPID && otherPID < np) {

      MPI_Send((void *)(&msgLength), 1, MPI_INT, otherPID, 0, MPI_COMM_WORLD);
      MPI_Send((void *)gMsg, msgLength, MPI_CHAR, otherPID, 1, MPI_COMM_WORLD);

    } else {
      opserr << "send -pid pid? data? - pid: " << otherPID << " invalid\n";
      return TCL_ERROR;
    }

  } else {
    if (myPID == 0) {
      MPI_Bcast((void *)(&msgLength), 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast((void *)gMsg, msgLength, MPI_CHAR, 0, MPI_COMM_WORLD);
    } else {
      opserr << "send data - only process 0 can do a broadcast - you may need "
                "to kill the application";
      return TCL_ERROR;
    }
  }
  return TCL_OK;
}

int
opsRecv(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
#ifdef _PARALLEL_MP
  MachineBroker* theMachineBroker = (MachineBroker*)clientData;
  if (argc < 2)
    return TCL_OK;

  int otherPID = 0;
  int myPID = theMachineBroker->getPID();
  int np = theMachineBroker->getNP();
  TCL_Char *varToSet = argv[argc - 1];

  int msgLength = 0;
  char *gMsg = 0;

  if (strcmp(argv[1], "-pid") == 0 && argc > 3) {

    bool fromAny = false;

    if ((strcmp(argv[2], "ANY") == 0) || (strcmp(argv[2], "ANY_SOURCE") == 0) ||
        (strcmp(argv[2], "MPI_ANY_SOURCE") == 0)) {
      fromAny = true;
    } else {
      if (Tcl_GetInt(interp, argv[2], &otherPID) != TCL_OK) {
        opserr << "recv -pid pid? data? - pid: " << argv[2] << " invalid\n";
        return TCL_ERROR;
      }
    }

    if (otherPID > -1 && otherPID < np) {
      MPI_Status status;

      if (fromAny == false)
        if (myPID != otherPID)
          MPI_Recv((void *)(&msgLength), 1, MPI_INT, otherPID, 0,
                   MPI_COMM_WORLD, &status);
        else {
          opserr << "recv -pid pid? data? - " << otherPID
                 << " cant receive from self!\n";
          return TCL_ERROR;
        }
      else {
        MPI_Recv((void *)(&msgLength), 1, MPI_INT, MPI_ANY_SOURCE, 0,
                 MPI_COMM_WORLD, &status);
        otherPID = status.MPI_SOURCE;
      }

      if (msgLength > 0) {
        gMsg = new char[msgLength];

        MPI_Recv((void *)gMsg, msgLength, MPI_CHAR, otherPID, 1, MPI_COMM_WORLD,
                 &status);

        Tcl_SetVar(interp, varToSet, gMsg, TCL_LEAVE_ERR_MSG);
      }

    } else {
      opserr << "recv -pid pid? data? - " << otherPID << " invalid\n";
      return TCL_ERROR;
    }
  } else {
    if (myPID != 0) {
      MPI_Bcast((void *)(&msgLength), 1, MPI_INT, 0, MPI_COMM_WORLD);

      if (msgLength > 0) {
        gMsg = new char[msgLength];

        MPI_Bcast((void *)gMsg, msgLength, MPI_CHAR, 0, MPI_COMM_WORLD);

        Tcl_SetVar(interp, varToSet, gMsg, TCL_LEAVE_ERR_MSG);
      }

    } else {
      opserr << "recv data - only process 0 can do a broadcast - you may need "
                "to kill the application";
      return TCL_ERROR;
    }
  }
#endif

  return TCL_OK;
}


int
opsPartition(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
#ifdef _PARALLEL_SP
  PartitionRuntime& part = *static_cast<PartitionRuntime*>(clientData);

  int eleTag;
  if (argc == 2) {
    if (Tcl_GetInt(interp, argv[1], &eleTag) != TCL_OK) {
      ;
    }
  }
  partitionModel(part, eleTag);
#endif
  return TCL_OK;
}

int 
wipePP(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

#ifdef _PARALLEL_SP

  PartitionRuntime& part = *static_cast<PartitionRuntime*>(clientData);

  if (OPS_PARTITIONED == true && part.num_subdomains > 1) {
    SubdomainIter &theSubdomains = part.theDomain.getSubdomains();
    Subdomain *theSub =nullptr;
    
    // create the appropriate domain decomposition analysis
    while ((theSub = theSubdomains()) != nullptr)
      theSub->wipeAnalysis();
  }
#endif
  return TCL_OK;  
}

