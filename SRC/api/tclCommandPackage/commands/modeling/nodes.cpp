/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
//
#include <assert.h>
#include <tcl.h>
#include <Node.h>
#include <Matrix.h>
#include <Domain.h>
#include <runtime/BasicModelBuilder.h>

int
TclCommand_addNode(ClientData clientData, Tcl_Interp *interp, int argc,
                   TCL_Char **argv)
{
  assert(clientData != nullptr);

  BasicModelBuilder *theTclBuilder = (BasicModelBuilder*)clientData;

  Domain *theTclDomain = theTclBuilder->getDomain();

  // ensure the destructor has not been called -
  if (theTclBuilder == 0 || clientData == 0) {
    opserr << "WARNING builder has been destroyed" << endln;
    return TCL_ERROR;
  }

  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  // make sure corect number of arguments on command line
  if (argc < 2 + ndm) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "        Want: node nodeTag? [ndm coordinates?] <-mass [ndf values?]>\n";
    return TCL_ERROR;
  }

  Node *theNode = 0;

  // get the nodal id
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << "WARNING invalid nodeTag\n";
    opserr << "Want: node nodeTag? [ndm coordinates?] <-mass [ndf values?]>\n";
    return TCL_ERROR;
  }

  // read in the coordinates and create the node
  double xLoc, yLoc, zLoc;
  if (ndm == 1) {
    // create a node in 1d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      opserr << "WARNING invalid XCoordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
  }

  else if (ndm == 2) {
    // create a node in 2d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      opserr << "WARNING invalid 1st coordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &yLoc) != TCL_OK) {
      opserr << "WARNING invalid 2nd coordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
  }

  else if (ndm == 3) {
    // create a node in 3d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      opserr << "WARNING invalid 1st coordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &yLoc) != TCL_OK) {
      opserr << "WARNING invalid 2nd coordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &zLoc) != TCL_OK) {
      opserr << "WARNING invalid 3rd coordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
  } else {
    opserr << "WARNING invalid ndm\n";
    opserr << "node: " << nodeId << endln;
    ;
    return TCL_ERROR;
  }

  // check for -ndf override option
  int currentArg = 2 + ndm;
  if (currentArg < argc && strcmp(argv[currentArg], "-ndf") == 0) {
    if (Tcl_GetInt(interp, argv[currentArg + 1], &ndf) != TCL_OK) {
      opserr << "WARNING invalid nodal ndf given for node " << nodeId << endln;
      return TCL_ERROR;
    }
    currentArg += 2;
  }

  //
  // create the node
  //
  if (ndm == 1)
    theNode = new Node(nodeId, ndf, xLoc);
  else if (ndm == 2)
    theNode = new Node(nodeId, ndf, xLoc, yLoc);
  else
    theNode = new Node(nodeId, ndf, xLoc, yLoc, zLoc);

  //
  // add the node to the domain
  //
  if (theTclDomain->addNode(theNode) == false) {
    opserr << "WARNING failed to add node to the domain\n";
    opserr << "node: " << nodeId << endln;
    delete theNode; // otherwise memory leak
    return TCL_ERROR;
  }

  while (currentArg < argc) {
    if (strcmp(argv[currentArg], "-mass") == 0) {
      currentArg++;
      if (argc < currentArg + ndf) {
        opserr << "WARNING incorrect number of nodal mass terms\n";
        opserr << "node: " << nodeId << endln;
        return TCL_ERROR;
      }
      Matrix mass(ndf, ndf);
      double theMass;
      for (int i = 0; i < ndf; i++) {
        if (Tcl_GetDouble(interp, argv[currentArg++], &theMass) != TCL_OK) {
          opserr << "WARNING invalid nodal mass term\n";
          opserr << "node: " << nodeId << ", dof: " << i + 1 << endln;
          return TCL_ERROR;
        }
        mass(i, i) = theMass;
      }
      theNode->setMass(mass);
    } else if (strcmp(argv[currentArg], "-dispLoc") == 0) {
      currentArg++;
      if (argc < currentArg + ndm) {
        opserr << "WARNING incorrect number of nodal display location terms, "
                  "need ndm\n";
        opserr << "node: " << nodeId << endln;
        return TCL_ERROR;
      }
      Vector displayLoc(ndm);
      double theCrd;
      for (int i = 0; i < ndm; i++) {
        if (Tcl_GetDouble(interp, argv[currentArg++], &theCrd) != TCL_OK) {
          opserr << "WARNING invalid nodal mass term\n";
          opserr << "node: " << nodeId << ", dof: " << i + 1 << endln;
          return TCL_ERROR;
        }
        displayLoc(i) = theCrd;
      }
      theNode->setDisplayCrds(displayLoc);

    } else if (strcmp(argv[currentArg], "-disp") == 0) {
      currentArg++;
      if (argc < currentArg + ndf) {
        opserr << "WARNING incorrect number of nodal disp terms\n";
        opserr << "node: " << nodeId << endln;
        return TCL_ERROR;
      }
      Vector disp(ndf);
      double theDisp;
      for (int i = 0; i < ndf; i++) {
        if (Tcl_GetDouble(interp, argv[currentArg++], &theDisp) != TCL_OK) {
          opserr << "WARNING invalid nodal disp term\n";
          opserr << "node: " << nodeId << ", dof: " << i + 1 << endln;
          return TCL_ERROR;
        }
        disp(i) = theDisp;
      }
      theNode->setTrialDisp(disp);
      theNode->commitState();

    } else if (strcmp(argv[currentArg], "-vel") == 0) {
      currentArg++;
      if (argc < currentArg + ndf) {
        opserr << "WARNING incorrect number of nodal vel terms\n";
        opserr << "node: " << nodeId << endln;
        return TCL_ERROR;
      }
      Vector disp(ndf);
      double theDisp;
      for (int i = 0; i < ndf; i++) {
        if (Tcl_GetDouble(interp, argv[currentArg++], &theDisp) != TCL_OK) {
          opserr << "WARNING invalid nodal vel term\n";
          opserr << "node: " << nodeId << ", dof: " << i + 1 << endln;
          return TCL_ERROR;
        }
        disp(i) = theDisp;
      }
      theNode->setTrialVel(disp);
      theNode->commitState();

    } else
      currentArg++;
  }

  // if we get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}

int
TclCommand_addNodalMass(ClientData clientData, Tcl_Interp *interp, int argc,
                        TCL_Char **argv)
{
  assert(clientData != nullptr);

  BasicModelBuilder *theTclBuilder = (BasicModelBuilder*)clientData;

  Domain *theTclDomain = theTclBuilder->getDomain();

  if (theTclBuilder == 0 || clientData == 0) {
    opserr << "WARNING builder has been destroyed - load \n";
    return TCL_ERROR;
  }

  int ndf = argc - 2;

  // make sure at least one other argument to contain type of system
  if (argc < (2 + ndf)) {
    opserr << "WARNING bad command - want: mass nodeId " << ndf << " mass values\n"; 
    return TCL_ERROR;
  }

  // get the id of the node
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << "WARNING invalid nodeId: " << argv[1];
    opserr << " - mass nodeId " << ndf << " forces\n";
    return TCL_ERROR;
  }

  // check for mass terms
  Matrix mass(ndf,ndf);
  double theMass;
  for (int i=0; i<ndf; i++)
  {
     if (Tcl_GetDouble(interp, argv[i+2], &theMass) != TCL_OK)
     {
          opserr << "WARNING invalid nodal mass term\n";
          opserr << "node: " << nodeId << ", dof: " << i+1 << endln;
          return TCL_ERROR;
      }
      mass(i,i) = theMass;
  }

  if (theTclDomain->setMass(mass, nodeId) != 0) {
    opserr << "WARNING failed to set mass at node " << nodeId << endln;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



int
TclCommand_getNDM(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  Domain *the_domain = builder->getDomain();

  int ndm;

  if (argc > 1) {
    int tag;
    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
      opserr << "WARNING ndm nodeTag? \n";
      return TCL_ERROR;
    }
    Node *theNode = the_domain->getNode(tag);
    if (theNode == 0) {
      opserr << "WARNING nodeTag " << tag << " does not exist \n";
      return TCL_ERROR;
    }
    const Vector &coords = theNode->getCrds();
    ndm = coords.Size();

  } else {
      ndm = builder->getNDM();
  }

  char buffer[20];
  sprintf(buffer, "%d", ndm);
  Tcl_AppendResult(interp, buffer, NULL);

  return TCL_OK;
}

int
TclCommand_getNDF(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  Domain *the_domain = builder->getDomain();
  int ndf;

  if (argc > 1) {
    int tag;
    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
      opserr << "WARNING ndf nodeTag? \n";
      return TCL_ERROR;
    }
    Node *theNode = the_domain->getNode(tag);
    if (theNode == 0) {
      opserr << "WARNING nodeTag " << tag << " does not exist \n";
      return TCL_ERROR;
    }
    ndf = theNode->getNumberDOF();

  } else {
      ndf = builder->getNDF();
  }

  char buffer[G3_NUM_DOF_BUFFER];
  if (abs(ndf) <  G3_MAX_NUM_DOFS){
    sprintf(buffer, "%d", ndf);
  } else {
    opserr << "ERROR -- Invalid DOF count encountered; got '" << ndf << "'.\n";
    return TCL_ERROR;
  }

  Tcl_AppendResult(interp, buffer, NULL);

  return TCL_OK;
}
