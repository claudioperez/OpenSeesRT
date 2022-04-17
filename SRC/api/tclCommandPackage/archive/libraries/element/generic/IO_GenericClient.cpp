

#include <g3_api.h>


#include <SRC/element/generic/GenericClient.h>
void *OPS_GenericClient()
{
  int ndf = OPS_GetNDF();
  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element genericClient eleTag -node Ndi Ndj ... -dof "
              "dofNdi -dof dofNdj ... -server ipPort <ipAddr> <-ssl> <-udp> "
              "<-dataSize size> <-noRayleigh>\n";
    return 0;
  }

  // tags
  int tag;
  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING: invalid tag\n";
    return 0;
  }

  // nodes
  const char *type = OPS_GetString();
  if (strcmp(type, "-node") != 0) {
    opserr << "WARNING expecting -node Ndi Ndj ...\n";
    return 0;
  }
  ID nodes(32);
  int numNodes = 0;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    int node;
    numdata = 1;
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (OPS_GetIntInput(&numdata, &node) < 0) {
      if (numArgs > OPS_GetNumRemainingInputArgs()) {
        // move current arg back by one
        OPS_ResetCurrentInputArg(-1);
      }
      break;
    }
    nodes(numNodes++) = node;
  }
  nodes.resize(numNodes);

  // dofs
  int numDOF = 0;
  ID *dofs = new ID[numNodes];
  for (int i = 0; i < numNodes; i++) {
    type = OPS_GetString();
    if (strcmp(type, "-dof") != 0 && strcmp(type, "-dir") != 0) {
      opserr << "WARNING expecting -dof dofNd" << i + 1 << ", but got " << type
             << endln;
      return 0;
    }
    ID dofsi(ndf);
    int numDOFi = 0;
    while (OPS_GetNumRemainingInputArgs() > 0) {
      int dof;
      numdata = 1;
      int numArgs = OPS_GetNumRemainingInputArgs();
      if (OPS_GetIntInput(&numdata, &dof) < 0) {
        if (numArgs > OPS_GetNumRemainingInputArgs()) {
          // move current arg back by one
          OPS_ResetCurrentInputArg(-1);
        }
        break;
      }
      if (dof < 1 || ndf < dof) {
        opserr << "WARNING invalid dof ID\n";
        return 0;
      }
      dofsi(numDOFi++) = dof - 1;
      numDOF++;
    }
    dofsi.resize(numDOFi);
    dofs[i] = dofsi;
  }

  // ipPort
  int ipPort;
  numdata = 1;
  type = OPS_GetString();
  if (strcmp(type, "-server") != 0) {
    opserr << "WARNING expecting -server ipPort <ipAddr>\n";
    return 0;
  }
  if (OPS_GetIntInput(&numdata, &ipPort) < 0) {
    opserr << "WARNING: invalid ipPort\n";
    return 0;
  }

  // options
  char *ipAddr = new char[10];
  strcpy(ipAddr, "127.0.0.1");
  int ssl = 0, udp = 0;
  int dataSize = 256;
  int doRayleigh = 1;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    type = OPS_GetString();
    if (strcmp(type, "-ssl") != 0 && strcmp(type, "-udp") != 0 &&
        strcmp(type, "-dataSize") != 0 && strcmp(type, "-noRayleigh") != 0 &&
        strcmp(type, "-doRayleigh") != 0) {
      delete[] ipAddr;
      ipAddr = new char[strlen(type) + 1];
      strcpy(ipAddr, type);
    } else if (strcmp(type, "-ssl") == 0) {
      ssl = 1;
      udp = 0;
    } else if (strcmp(type, "-udp") == 0) {
      udp = 1;
      ssl = 0;
    } else if (strcmp(type, "-dataSize") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING wrong dataSize specified\n";
        return 0;
      }
      numdata = 1;
      if (OPS_GetIntInput(&numdata, &dataSize) < 0) {
        opserr << "WARNING invalid dataSize value\n";
        return 0;
      }
    } else if (strcmp(type, "-doRayleigh") == 0) {
      doRayleigh = 1;
    } else if (strcmp(type, "-noRayleigh") == 0) {
      doRayleigh = 0;
    }
  }

  // create object
  Element *theEle = new GenericClient(tag, nodes, dofs, ipPort, ipAddr, ssl,
                                      udp, dataSize, doRayleigh);

  // cleanup dynamic memory
  if (dofs != 0)
    delete[] dofs;
  delete[] ipAddr;

  return theEle;
}
