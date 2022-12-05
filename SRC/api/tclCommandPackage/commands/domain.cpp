#include <assert.h>
#include <set>
#include <vector>
#include <algorithm>     // std::sort

#include <g3_api.h>
#include <FileStream.h>
#include <SimulationInformation.h>

#include <Domain.h>
#include <LoadPattern.h>
#include <Parameter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>

#include <Pressure_Constraint.h>
#include <Element.h>
#include <ElementIter.h>

#include <Node.h>
#include <NodeIter.h>
#include <Vector.h>

extern  SimulationInformation simulationInfo;

int printElement(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char **argv, OPS_Stream &output);

int printNode(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char **argv, OPS_Stream &output);

int printIntegrator(ClientData clientData, Tcl_Interp *interp, int argc,
                    TCL_Char **argv, OPS_Stream &output);

int printAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc,
                   TCL_Char **argv, OPS_Stream &output);

int
printModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain * domain = (Domain*)clientData;
  // G3_Runtime *rt = G3_getRuntime(interp);
  // Domain *domain = G3_getDomain(rt);

  int currentArg = 1;
  int res = 0;

  int flag = OPS_PRINT_CURRENTSTATE;

  FileStream outputFile;
  OPS_Stream *output = &opserr;
  bool done = false;

  // if just 'print' then print out the entire domain
  if (argc == currentArg) {
    opserr << domain;
    return TCL_OK;
  }

  while (done == false) {
    // if 'print ele i j k..' print out some elements
    if ((strcmp(argv[currentArg], "-ele") == 0) ||
        (strcmp(argv[currentArg], "ele") == 0)) {
      currentArg++;
      res = printElement(clientData, interp, argc - currentArg, argv + currentArg, *output);
      done = true;
    }
    // if 'print node i j k ..' print out some nodes
    else if ((strcmp(argv[currentArg], "-node") == 0) ||
             (strcmp(argv[currentArg], "node") == 0)) {
      currentArg++;
      res = printNode(clientData, interp, argc - currentArg, argv + currentArg,
                      *output);
      done = true;
    }

    // if 'print integrator flag' print out the integrator
    else if ((strcmp(argv[currentArg], "integrator") == 0) ||
             (strcmp(argv[currentArg], "-integrator") == 0)) {
      currentArg++;
      res = printIntegrator(clientData, interp, argc - currentArg,
                            argv + currentArg, *output);
      done = true;
    }

    // if 'print algorithm flag' print out the algorithm
    else if ((strcmp(argv[currentArg], "algorithm") == 0) ||
             (strcmp(argv[currentArg], "-algorithm") == 0)) {
      currentArg++;
      res = printAlgorithm(clientData, interp, argc - currentArg,
                           argv + currentArg, *output);
      done = true;
    }
    else if ((strcmp(argv[currentArg], "-JSON") == 0)) {
      currentArg++;
      flag = OPS_PRINT_PRINTMODEL_JSON;
    }

    else {
      if ((strcmp(argv[currentArg], "file") == 0) ||
          (strcmp(argv[currentArg], "-file") == 0))
        currentArg++;

      openMode mode = APPEND;
      if (flag == OPS_PRINT_PRINTMODEL_JSON)
        mode = OVERWRITE;
      if (outputFile.setFile(argv[currentArg], mode) != 0) {
        opserr << "print <filename> .. - failed to open file: "
               << argv[currentArg] << endln;
        return TCL_ERROR;
      }
      currentArg++;

      // if just 'print <filename>' then print out the entire domain to eof
      if (argc == currentArg) {
        if (flag == OPS_PRINT_PRINTMODEL_JSON)
          simulationInfo.Print(outputFile, flag);

        domain->Print(outputFile, flag);
        return TCL_OK;
      }

      output = &outputFile;
    }
  }

  // close the output file
  outputFile.close();
  return res;
}

int
printElement(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv, OPS_Stream &output)
{
  assert(clientData != nullptr);
  Domain * the_domain = (Domain*)clientData;
  // Domain *the_domain = G3_getDomain(G3_getRuntime(interp));

  int flag = 0; // default flag sent to a nodes Print() method
  int eleArg = 0;

  // if just 'print <filename> node' print all the nodes - no flag
  if (argc == 0) {
    ElementIter &theElements = the_domain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0)
      theElement->Print(output);
    return TCL_OK;
  }

  // if 'print <filename> Element flag int <int int ..>' get the flag
  if ((strcmp(argv[0], "flag") == 0) ||
      (strcmp(argv[0], "-flag")) == 0) { // get the specified flag
    if (argc < 2) {
      opserr << "WARNING print <filename> ele <flag int> no int specified \n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[1], &flag) != TCL_OK) {
      opserr << "WARNING print ele failed to get integer flag: \n";
      opserr << argv[eleArg] << endln;
      return TCL_ERROR;
    }
    eleArg += 2;
  }

  // now print the Elements with the specified flag, 0 by default
  if (argc == eleArg) {
    ElementIter &theElements = the_domain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0)
      theElement->Print(output, flag);
    return TCL_OK;
  } else {

    // otherwise print out the specified nodes i j k .. with flag
    int numEle = argc - eleArg;
    ID *theEle = new ID(numEle);
    for (int i = 0; i < numEle; i++) {
      int eleTag;
      if (Tcl_GetInt(interp, argv[i + eleArg], &eleTag) != TCL_OK) {
        opserr << "WARNING print ele failed to get integer: " << argv[i]
               << endln;
        return TCL_ERROR;
      }
      (*theEle)(i) = eleTag;
    }

    the_domain->Print(output, 0, theEle, flag);
    delete theEle;
  }

  return TCL_OK;
}

// function to print out the nodal information conatined in line
//     print <filename> node <flag int> <int int int>
// Parameters
//   nodeArg: integer equal to arg count to node plus 1
//   output:  output stream to which the results are sent
//
int
printNode(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv,
          OPS_Stream &output)
{
  int flag = 0; // default flag sent to a nodes Print() method
  int nodeArg = 0;
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain* domain = G3_getDomain(rt);

  // if just 'print <filename> node' print all the nodes - no flag
  if (argc == 0) {
    NodeIter &theNodes = domain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0)
      theNode->Print(output);
    return TCL_OK;
  }

  // if 'print <filename> node flag int <int int ..>' get the flag
  if ((strcmp(argv[0], "flag") == 0) || (strcmp(argv[0], "-flag") == 0)) {
    // get the specified flag
    if (argc <= nodeArg) {
      opserr << "WARNING print <filename> node <flag int> no int specified \n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[1], &flag) != TCL_OK) {
      opserr << "WARNING print node failed to get integer flag: \n";
      opserr << argv[nodeArg] << endln;
      return TCL_ERROR;
    }
    nodeArg += 2;
  }

  // now print the nodes with the specified flag, 0 by default

  // if 'print <filename> node flag'
  //     print out all the nodes in the domain with flag
  if (nodeArg == argc) {
    NodeIter &theNodes = domain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0)
      theNode->Print(output, flag);
    return TCL_OK;
  } else {
    // otherwise print out the specified nodes i j k .. with flag
    int numNodes = argc - nodeArg;
    ID *theNodes = new ID(numNodes);
    for (int i = 0; i < numNodes; i++) {
      int nodeTag;
      if (Tcl_GetInt(interp, argv[nodeArg], &nodeTag) != TCL_OK) {
        opserr << "WARNING print node failed to get integer: " << argv[nodeArg]
               << endln;
        return TCL_ERROR;
      }
      (*theNodes)(i) = nodeTag;
      nodeArg++;
    }

    domain->Print(output, theNodes, 0, flag);
    delete theNodes;
  }

  return TCL_OK;
}

int
removeObject(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain * the_domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - remove objectType?\n";
    return TCL_ERROR;
  }

  int tag;
  if ((strcmp(argv[1], "element") == 0) || (strcmp(argv[1], "ele") == 0)) {
    if (argc < 3) {
      opserr << "WARNING want - remove element eleTag?\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING remove element tag? failed to read tag: " << argv[2]
             << endln;
      return TCL_ERROR;
    }
    Element *theEle = the_domain->removeElement(tag);
    if (theEle != 0) {
#if 0
      // we also have to remove any elemental loads from the domain
      LoadPatternIter &theLoadPatterns = the_domain->getLoadPatterns();
      LoadPattern *thePattern;

      // go through all load patterns
      while ((thePattern = theLoadPatterns()) != 0) {
        ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
        ElementalLoad *theLoad;

        // go through all elemental loads in the pattern
        while ((theLoad = theEleLoads()) != 0) {

          // remove & destroy elemental from elemental load if there
          // note - if last element in load, remove the load and delete it

          /* *****************
             int numLoadsLeft = theLoad->removeElement(tag);
             if (numLoadsLeft == 0) {
             thePattern->removeElementalLoad(theLoad->getTag());
             delete theLoad;
             }
          *********************/
        }
      }
#endif
      // finally invoke the destructor on the element
      delete theEle;
    }
  }

  else if (strcmp(argv[1], "loadPattern") == 0) {
    if (argc < 3) {
      opserr << "WARNING want - remove loadPattern patternTag?\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING remove loadPattern tag? failed to read tag: "
             << argv[2] << endln;
      return TCL_ERROR;
    }
    LoadPattern *thePattern = the_domain->removeLoadPattern(tag);
    if (thePattern != 0) {
      thePattern->clearAll();
      delete thePattern;
    }
  }
#if 0
  else if ((strcmp(argv[1], "TimeSeries") == 0) ||
           (strcmp(argv[1], "timeSeries") == 0)) {
    if (argc < 3) {
      opserr << "WARNING want - remove loadPattern patternTag?\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING remove loadPattern tag? failed to read tag: "
             << argv[2] << endln;
      return TCL_ERROR;
    }
    bool ok = OPS_removeTimeSeries(tag);
    if (ok == true)
      return TCL_OK;
    else
      return TCL_ERROR;
  }
#endif
  else if (strcmp(argv[1], "parameter") == 0) {
    if (argc < 3) {
      opserr << "WARNING want - remove parameter paramTag?\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING remove parameter tag? failed to read tag: " << argv[2]
             << endln;
      return TCL_ERROR;
    }
    Parameter *theParameter = the_domain->removeParameter(tag);
    if (theParameter != 0) {
      delete theParameter;
    }
  }

  else if (strcmp(argv[1], "node") == 0) {
    if (argc < 3) {
      opserr << "WARNING want - remove node nodeTag?\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING remove node tag? failed to read tag: " << argv[2]
             << endln;
      return TCL_ERROR;
    }
    Node *theNode = the_domain->removeNode(tag);
    if (theNode != 0) {
      delete theNode;
    }
    Pressure_Constraint *thePC = the_domain->removePressure_Constraint(tag);
    if (thePC != 0) {
      delete thePC;
    }
  }

  else if (strcmp(argv[1], "recorders") == 0) {
    the_domain->removeRecorders();
  }

  else if ((strcmp(argv[1], "recorder") == 0)) {
    if (argc < 3) {
      opserr << "WARNING want - remove recorder recorderTag?\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING remove recorder tag? failed to read tag: " << argv[2]
             << endln;
      return TCL_ERROR;
    }
    return the_domain->removeRecorder(tag);
  }


  else if ((strcmp(argv[1], "SPconstraint") == 0) ||
           (strcmp(argv[1], "sp") == 0)) {
    if (argc < 3) {
      opserr << "WARNING want - remove SPconstraint spTag? -or- remove "
                "SPconstraint nodeTag? dofTag? <patternTag?>\n";
      return TCL_ERROR;
    }
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
        opserr << "WARNING remove sp tag? failed to read tag: " << argv[2]
               << endln;
        return TCL_ERROR;
      }
      SP_Constraint *theSPconstraint = the_domain->removeSP_Constraint(tag);
      if (theSPconstraint != 0) {
        delete theSPconstraint;
      }
    } else {
      int nodeTag, dofTag;
      int patternTag = -1;

      if (Tcl_GetInt(interp, argv[2], &nodeTag) != TCL_OK) {
        opserr << "WARNING remove sp tag? failed to read node tag: " << argv[2]
               << endln;
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[3], &dofTag) != TCL_OK) {
        opserr << "WARNING remove sp tag? failed to read dof tag: " << argv[3]
               << endln;
        return TCL_ERROR;
      }

      if (argc == 5) {
        if (Tcl_GetInt(interp, argv[4], &patternTag) != TCL_OK) {
          opserr << "WARNING remove sp tag? failed to read pattern tag: "
                 << argv[4] << endln;
          return TCL_ERROR;
        }
      }
      dofTag--; // one for C++ indexing of dof

      the_domain->removeSP_Constraint(nodeTag, dofTag, patternTag);

      return TCL_OK;
    }
  }

  else if ((strcmp(argv[1], "MPconstraint") == 0) ||
           (strcmp(argv[1], "mp") == 0)) {
    if (argc < 3) {
      opserr << "WARNING want - remove MPconstraint nNodeTag? -or- remove "
                "MPconstraint -tag mpTag\n";
      return TCL_ERROR;
    }
    int nodTag = 0;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &nodTag) != TCL_OK) {
        opserr << "WARNING remove mp nodeTag? failed to read nodeTag: "
               << argv[2] << endln;
        return TCL_ERROR;
      }

      the_domain->removeMP_Constraints(nodTag);
      return TCL_OK;
    }
    if (strcmp(argv[2], "-tag") == 0 && argc > 3) {
      if (Tcl_GetInt(interp, argv[3], &nodTag) != TCL_OK) {
        opserr << "WARNING remove mp -tag mpTag? failed to read mpTag: "
               << argv[3] << endln;
        return TCL_ERROR;
      }

      the_domain->removeMP_Constraint(nodTag);
      return TCL_OK;
    }
  }

#ifdef _RELIABILITY
  // AddingSensitivity:BEGIN ///////////////////////////////////////
  else if (strcmp(argv[1], "randomVariable") == 0) {
    int rvTag;
    if (Tcl_GetInt(interp, argv[2], &rvTag) != TCL_OK) {
      opserr << "WARNING invalid input: rvTag \n";
      return TCL_ERROR;
    }
    ReliabilityDomain *theReliabilityDomain =
        theReliabilityBuilder->getReliabilityDomain();
    theReliabilityDomain->removeRandomVariable(rvTag);
  } else if (strcmp(argv[1], "performanceFunction") == 0) {
    int lsfTag;
    if (Tcl_GetInt(interp, argv[2], &lsfTag) != TCL_OK) {
      opserr << "WARNING invalid input: lsfTag \n";
      return TCL_ERROR;
    }
    ReliabilityDomain *theReliabilityDomain =
        theReliabilityBuilder->getReliabilityDomain();
    theReliabilityDomain->removeLimitStateFunction(lsfTag);
  } else if (strcmp(argv[1], "cutset") == 0) {
    int cutTag;
    if (Tcl_GetInt(interp, argv[2], &cutTag) != TCL_OK) {
      opserr << "WARNING invalid input: cutTag \n";
      return TCL_ERROR;
    }
    ReliabilityDomain *theReliabilityDomain =
        theReliabilityBuilder->getReliabilityDomain();
    theReliabilityDomain->removeCutset(cutTag);
  } else if (strcmp(argv[1], "sensitivityAlgorithm") == 0) {
    if (theSensitivityAlgorithm != 0) {
      // the_static_analysis->setSensitivityAlgorithm(0);
      theSensitivityAlgorithm = 0;
      theSensitivityIntegrator = 0;
    }
  }
// AddingSensitivity:END ///////////////////////////////////////
#endif

  else
    opserr << "WARNING remove " << argv[1] << " not supported" << endln;

  return TCL_OK;
}


int
fixedNodes(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain * domain = (Domain*)clientData;

  SP_Constraint *theSP;
  SP_ConstraintIter &spIter = domain->getDomainAndLoadPatternSPs();

  // get unique constrained nodes with set
  std::set<int> tags;
  int tag;
  while ((theSP = spIter()) != 0) {
    tag = theSP->getNodeTag();
    tags.insert(tag);
  }
  // assign set to vector and sort
  std::vector<int> tagv;
  tagv.assign(tags.begin(), tags.end());
  std::sort(tagv.begin(), tagv.end());
  // loop through unique, sorted tags, adding to output
  char buffer[20];
  for (int tag : tagv) {
    sprintf(buffer, "%d ", tag);
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
fixedDOFs(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain * theDomain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING want - fixedDOFs fNode?\n";
    return TCL_ERROR;
  }

  int fNode;
  if (Tcl_GetInt(interp, argv[1], &fNode) != TCL_OK) {
    opserr << "WARNING fixedDOFs fNode? - could not read fNode? \n";
    return TCL_ERROR;
  }

  SP_Constraint *theSP;
  SP_ConstraintIter &spIter = theDomain->getDomainAndLoadPatternSPs();

  int tag;
  Vector fixed(6);
  while ((theSP = spIter()) != 0) {
    tag = theSP->getNodeTag();
    if (tag == fNode) {
      fixed(theSP->getDOF_Number()) = 1;
    }
  }

  char buffer[20];
  for (int i = 0; i < 6; i++) {
    if (fixed(i) == 1) {
      sprintf(buffer, "%d ", i + 1);
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}

int
constrainedNodes(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain * theDomain = (Domain*)clientData;

  bool all = 1;
  int rNode;
  if (argc > 1) {
    if (Tcl_GetInt(interp, argv[1], &rNode) != TCL_OK) {
      opserr << "WARNING constrainedNodes <rNode?> - could not read rNode? \n";
      return TCL_ERROR;
    }
    all = 0;
  }

  MP_Constraint *theMP;
  MP_ConstraintIter &mpIter = theDomain->getMPs();

  // get unique constrained nodes with set
  std::set<int> tags;
  int tag;
  while ((theMP = mpIter()) != 0) {
    tag = theMP->getNodeConstrained();
    if (all || rNode == theMP->getNodeRetained()) {
      tags.insert(tag);
    }
  }
  // assign set to vector and sort
  std::vector<int> tagv;
  tagv.assign(tags.begin(), tags.end());
  std::sort(tagv.begin(), tagv.end());
  // loop through unique, sorted tags, adding to output
  char buffer[20];
  for (int tag : tagv) {
    sprintf(buffer, "%d ", tag);
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
constrainedDOFs(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain * theDomain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING want - constrainedDOFs cNode? <rNode?> <rDOF?>\n";
    return TCL_ERROR;
  }

  int cNode;
  if (Tcl_GetInt(interp, argv[1], &cNode) != TCL_OK) {
    opserr << "WARNING constrainedDOFs cNode? <rNode?> <rDOF?> - could not "
              "read cNode? \n";
    return TCL_ERROR;
  }

  int rNode;
  bool allNodes = 1;
  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &rNode) != TCL_OK) {
      opserr << "WARNING constrainedDOFs cNode? <rNode?> <rDOF?> - could not "
                "read rNode? \n";
      return TCL_ERROR;
    }
    allNodes = 0;
  }

  int rDOF;
  bool allDOFs = 1;
  if (argc > 3) {
    if (Tcl_GetInt(interp, argv[3], &rDOF) != TCL_OK) {
      opserr << "WARNING constrainedDOFs cNode? <rNode?> <rDOF?> - could not "
                "read rDOF? \n";
      return TCL_ERROR;
    }
    rDOF--;
    allDOFs = 0;
  }

  MP_Constraint *theMP;
  MP_ConstraintIter &mpIter = theDomain->getMPs();

  int tag;
  int i;
  int n;
  Vector constrained(6);
  while ((theMP = mpIter()) != 0) {
    tag = theMP->getNodeConstrained();
    if (tag == cNode) {
      if (allNodes || rNode == theMP->getNodeRetained()) {
        const ID &cDOFs = theMP->getConstrainedDOFs();
        n = cDOFs.Size();
        if (allDOFs) {
          for (i = 0; i < n; i++) {
            constrained(cDOFs(i)) = 1;
          }
        } else {
          const ID &rDOFs = theMP->getRetainedDOFs();
          for (i = 0; i < n; i++) {
            if (rDOF == rDOFs(i))
              constrained(cDOFs(i)) = 1;
          }
        }
      }
    }
  }
  char buffer[20];
  for (int i = 0; i < 6; i++) {
    if (constrained(i) == 1) {
      sprintf(buffer, "%d ", i + 1);
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}

