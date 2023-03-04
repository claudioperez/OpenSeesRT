/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Description: Domain manipulation commands which do not require
// additional namespacing.
//
#include <assert.h>
#include <set>
#include <vector>
#include <algorithm>     // std::sort

#include <tcl.h>
#include <FileStream.h>

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


int
domainChange(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  ((Domain*)clientData)->domainChange();
  return TCL_OK;
}


int
removeObject(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain * the_domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of object
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
fixedNodes(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
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
fixedDOFs(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
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
                 TCL_Char ** const argv)
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
                TCL_Char ** const argv)
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

