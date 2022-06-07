#include <Block2D.h>
#include <Block3D.h>

int
TclCommand_doBlock2D(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char **argv)
{ 
  G3_Runtime* rt = G3_getRuntime(interp);
  Domain *theTclDomain = G3_getDomain(rt);
  int ndm = G3_getNDM(rt);
  if (ndm < 2) {
    opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    opserr << " : model dimension (ndm) must be at leat 2 " << endln;
    return TCL_ERROR;
  }

  if (argc < 8) {
    opserr << "WARNING incorrect numer of args :block2D numX? numY? startNode? startEle? eleType? eleArgs?"; 
    return TCL_ERROR;
  }
  int numX, numY, startNodeNum, startEleNum;
  if (Tcl_GetInt(interp, argv[1], &numX) != TCL_OK) {
    opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid numX: " << argv[1] << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &numY) != TCL_OK) {
    opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid numY: " << argv[2] << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &startNodeNum) != TCL_OK) {
    opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid startNode: " << argv[3] << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[4], &startEleNum) != TCL_OK) {
    opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid startEle: " << argv[4] << endln;
    return TCL_ERROR;
  }


  static Matrix Coordinates(9,3);
  static ID     haveNode(9);
  Coordinates.Zero();
  for (int k=0; k<9; k++) haveNode(k) = -1;

  int numNodes = 4;
  if (argc == 10) {
    if (strcmp(argv[7],"-numEleNodes") == 0)
      if (Tcl_GetInt(interp, argv[8], &numNodes) != TCL_OK) {
        opserr << "WARNING block2D numX? numY? startNode? startEle? eleType eleArgs?";
        opserr << " -numEleNodes numNodes?: invalid numNodes: " << argv[8] << endln; 
        return TCL_ERROR;
      }
    if (numNodes != 4 && numNodes != 9) {
      opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs? ";
      opserr << "-numEleNodes numNodes?: invalid numNodes: " << argv[8] << " 4 or 9 only\n";
      return TCL_ERROR;
    }

    if (numNodes == 9) {
      if (((numX % 2) != 0) || ((numY % 2) != 0)) {
        opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs? ";
        opserr << "-numEleNodes 9: numX and numY MUST BOTH BE EVEN\n";
        return TCL_ERROR;
      }
    }
  }

  int nodalInfo = 9;
  if (numNodes == 4)
    nodalInfo = 7;

  TCL_Char **argvNodes;
  int  argcNodes;

  Tcl_SplitList(interp, argv[nodalInfo], &argcNodes, &argvNodes);

  int ndf = G3_getNDF(rt);

  int count = 0;
  while (count < argcNodes) {
    int nodeTag;
    double value;
    if ((count + ndm + 1) >  argcNodes) {
      opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
      opserr << " : invalid number of node args: " << argv[7] << endln;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argvNodes[count], &nodeTag) != TCL_OK) {
      opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
      opserr << " : invalid node tag: " << argvNodes[count] << endln;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR;
    }
    if (nodeTag < 1 || nodeTag > 9) {
      opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
      opserr << " : invalid node tag out of bounds [1,9]: " << argvNodes[count] << endln;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR;
    }
    for (int i=0; i<ndm; i++) {
      if (Tcl_GetDouble(interp, argvNodes[count+1+i], &value) != TCL_OK) {
        opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
        opserr << " : invalid node coordinate for node: " << argvNodes[count] << endln;
        Tcl_Free((char *)argvNodes); return TCL_ERROR;
      }
      Coordinates(nodeTag-1,i) = value;
      haveNode(nodeTag-1) = nodeTag;
    }
    count += 1 + ndm;
  }

  Tcl_Free((char *)argvNodes);

  Block2D  theBlock(numX, numY, haveNode, Coordinates, numNodes);

  // create the nodes: (numX+1)*(numY+1) nodes to be created
  int nodeID = startNodeNum;
  int jj;
  for (jj=0; jj<=numY; jj++) {
    for (int ii=0; ii<=numX; ii++) {
      const Vector &nodeCoords = theBlock.getNodalCoords(ii,jj);
      double xLoc = nodeCoords(0);
      double yLoc = nodeCoords(1);
      double zLoc = nodeCoords(2);
      Node *theNode = 0;
      if (ndm == 2) {
        theNode = new Node(nodeID,ndf,xLoc, yLoc);
      } else if (ndm == 3) {
        theNode = new Node(nodeID,ndf,xLoc, yLoc, zLoc);
      }
      if (theNode == 0) {
        opserr << "WARNING ran out of memory creating node\n";
        opserr << "node: " << nodeID << endln;
        return TCL_ERROR;
      }
      if (theTclDomain->addNode(theNode) == false) {
        opserr << "WARNING failed to add node to the domain\n";
        opserr << "node: " << nodeID << endln;
        delete theNode; // otherwise memory leak
        return TCL_ERROR;
      }
      nodeID++;
    }
  }

  // create the elements: numX*numY elements to be created if 4 node elements
  //                      numX/2 * numY /2 nodes to be v=created if 9 node elements 
  TCL_Char *eleType = argv[5]; 
  TCL_Char *additionalEleArgs = argv[6];
  //  const ID &nodeTags = theBlock.getElementNodes(0,0);
  //  int numNodes = nodeTags.Size();

  // assumes 15 is largest string for individual nodeTags
  count = int(10 + strlen(eleType) + strlen(additionalEleArgs) + 15 *(numNodes+1));
  char *eleCommand = new char[count]; int initialCount = int(8 + strlen(eleType));

  int  eleID = startEleNum;
  if (numNodes == 9) {
    numX /= 2;
    numY /= 2;
  }

  for (jj=0; jj<numY; jj++) {
    for (int ii=0; ii<numX; ii++) {
      count = initialCount;
      const ID &nodeTags = theBlock.getElementNodes(ii,jj);
      // create the string to be evaluated
      strcpy(eleCommand, "element ");
      strcpy(&eleCommand[8], eleType);
      count += sprintf(&eleCommand[count], " %d ", eleID);
      for (int i=0; i<numNodes; i++) {
        int nodeTag = nodeTags(i)+startNodeNum;
        count += sprintf(&eleCommand[count], " %d ", nodeTag);
      }
      strcat(eleCommand, additionalEleArgs);

      // now to create the element we get the string eveluated
      if (Tcl_Eval(interp, eleCommand) != TCL_OK) {
          delete [] eleCommand;
        return TCL_ERROR;
      }
      eleID++;
    }
  }
  delete [] eleCommand;
  return TCL_OK;
}


int
TclCommand_doBlock3D(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char **argv)
{
  G3_Runtime* rt = G3_getRuntime(interp);
  int ndm = G3_getNDM(rt);
  Domain *theTclDomain = G3_getDomain(rt);

  if (ndm < 3) {
    opserr << "WARNING block3D numX? numY? startNode? startEle? eleType? eleArgs?";
    opserr << " : model dimension (ndm) must be at leat 2 " << endln;
    return TCL_ERROR;
  }

  int numX, numY, numZ, startNodeNum, startEleNum;
  if (Tcl_GetInt(interp, argv[1], &numX) != TCL_OK) {
    opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid numX: " << argv[1] << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &numY) != TCL_OK) {
    opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid numY: " << argv[2] << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &numZ) != TCL_OK) {
    opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid numZ: " << argv[3] << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[4], &startNodeNum) != TCL_OK) {
    opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid startNode: " << argv[4] << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[5], &startEleNum) != TCL_OK) {
    opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid startEle: " << argv[5] << endln;
    return TCL_ERROR;
  }

  static Matrix Coordinates(27,3);
  static ID     haveNode(27);
  Coordinates.Zero();
  for (int k=0; k<27; k++) haveNode(k) = -1;

  TCL_Char *nodalInfo = argv[8];
  TCL_Char **argvNodes;
  int  argcNodes;

  Tcl_SplitList(interp, nodalInfo, &argcNodes, &argvNodes);

  int ndf = G3_getNDF(rt);

  int count = 0;
  while (count < argcNodes) {
    if ((count + ndm + 1) > argcNodes) {
      opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType eleArgs?";
      opserr << " : invalid number of node args: " << argv[8] << endln;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR;
    }
    int nodeTag;
    double value;
    if (Tcl_GetInt(interp, argvNodes[count], &nodeTag) != TCL_OK) {
      opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
      opserr << " : invalid node id in node args: " << argvNodes[count] << endln;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR;
    }
    if (nodeTag < 1 || nodeTag > 27) {
      opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
      opserr << " : node tag out of bounds [1, 27]: " << argvNodes[count] << endln;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR;
    }
    for (int i=0; i<ndm; i++) {
      if (Tcl_GetDouble(interp, argvNodes[count+1+i], &value) != TCL_OK) {
        opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
        opserr << " : invalid coordinate in node args: " << argvNodes[count] << endln;
        Tcl_Free((char *)argvNodes);
        return TCL_ERROR;
      }
      Coordinates(nodeTag-1,i) = value;
      haveNode(nodeTag-1) = nodeTag;
    }
    count += 1 + ndm;
  }

  Tcl_Free((char *)argvNodes);

  Block3D  theBlock(numX, numY, numZ, haveNode, Coordinates);

  // create the nodes: (numX+1)*(numY+1) nodes to be created
  int nodeID = startNodeNum;
  int kk;
  for (kk=0; kk<=numZ; kk++) {
    for (int jj=0; jj<=numY; jj++) {
      for (int ii=0; ii<=numX; ii++) {
        const Vector &nodeCoords = theBlock.getNodalCoords(ii,jj,kk);
        double xLoc = nodeCoords(0);
        double yLoc = nodeCoords(1);
        double zLoc = nodeCoords(2);
        Node *theNode = 0;
        theNode = new Node(nodeID,ndf,xLoc, yLoc, zLoc);

        if (theNode == 0) {
          opserr << "WARNING ran out of memory creating node\n";
          opserr << "node: " << nodeID << endln;
          return TCL_ERROR;
        }

        if (theTclDomain->addNode(theNode) == false) {
          opserr << "WARNING failed to add node to the domain\n";
          opserr << "node: " << nodeID << endln;
          delete theNode; // otherwise memory leak
          return TCL_ERROR;
        }

        nodeID++;
      }
    }
  }

  // create the elements: numX*numY elements to be created
  TCL_Char *eleType = argv[6];
  TCL_Char *additionalEleArgs = argv[7];
  const ID &nodeTags = theBlock.getElementNodes(0,0,0);
  int numNodes = nodeTags.Size();

  // assumes 15 is largest string for individual nodeTags
  count = int(10 + strlen(eleType) + strlen(additionalEleArgs) + 15 * (numNodes+1));
  char *eleCommand = new char[count];
  int initialCount = int(8 + strlen(eleType));

  int  eleID = startEleNum;
  for (kk=0; kk<numZ; kk++) {
    for (int jj=0; jj<numY; jj++) {
      for (int ii=0; ii<numX; ii++) {
        count = initialCount;

        const ID &nodeTags = theBlock.getElementNodes(ii,jj,kk);

        // create the string to be evaluated
        strcpy(eleCommand, "element ");
        strcpy(&eleCommand[8], eleType);
        count += sprintf(&eleCommand[count], " %d ", eleID);
        for (int i=0; i<numNodes; i++) {
          int nodeTag = nodeTags(i)+startNodeNum;
          count += sprintf(&eleCommand[count], " %d ", nodeTag);
        }
        strcat(eleCommand, additionalEleArgs);

        // now to create the element we get the string eveluated
        if (Tcl_Eval(interp, eleCommand) != TCL_OK) {
        delete [] eleCommand;
          return TCL_ERROR;
        }
        eleID++;
      }
    }
  }

  delete [] eleCommand;
  return TCL_OK;
}
