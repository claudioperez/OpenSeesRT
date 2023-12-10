/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Written: Minjie Zhu, Claudio Perez
//
#include "BasicAnalysisBuilder.h"
#include <Domain.h>
#include <assert.h>
#include <stdio.h>

#include <G3_Logging.h>

#include <EquiSolnAlgo.h>
#include <IncrementalIntegrator.h>
#include <StaticIntegrator.h>
#include <TransientIntegrator.h>
#include <LinearSOE.h>
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>
#include <DOF_Numberer.h>
#include <ConstraintHandler.h>
#include <ConvergenceTest.h>
#include <AnalysisModel.h>
#include <TimeSeries.h>
#include <LoadPattern.h>

// For eigen()
#include <FE_EleIter.h>
#include <FE_Element.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>

// Defaults
#include <Newmark.h>
#include <EigenSOE.h>
#include <SymBandEigenSolver.h>
#include <SymBandEigenSOE.h>
#include <FullGenEigenSolver.h>
#include <FullGenEigenSOE.h>
#include <ArpackSOE.h>
#include <ProfileSPDLinSOE.h>
#include <NewtonRaphson.h>
#include <RCM.h>
#include <LoadControl.h>
#include <ProfileSPDLinSolver.h>
#include <CTestNormUnbalance.h>
#include <ProfileSPDLinDirectSolver.h>
#include <PlainHandler.h>
#include <TransformationConstraintHandler.h>


BasicAnalysisBuilder::BasicAnalysisBuilder()
:theHandler(nullptr),theNumberer(nullptr),theAlgorithm(nullptr),
 theSOE(nullptr),theEigenSOE(nullptr),theStaticIntegrator(nullptr),theTransientIntegrator(nullptr),
 theTest(nullptr),theStaticAnalysis(nullptr),theTransientAnalysis(nullptr),
 theVariableTimeStepTransientAnalysis(nullptr),CurrentAnalysisFlag(CURRENT_EMPTY_ANALYSIS)
{
  theAnalysisModel = new AnalysisModel();
}

BasicAnalysisBuilder::BasicAnalysisBuilder(Domain* domain)
:theHandler(nullptr),theNumberer(nullptr),theAlgorithm(nullptr),
 theSOE(nullptr),theEigenSOE(nullptr),theStaticIntegrator(nullptr),theTransientIntegrator(nullptr),
 theTest(nullptr),theStaticAnalysis(nullptr),theTransientAnalysis(nullptr),
 theVariableTimeStepTransientAnalysis(nullptr), theDomain(domain),
 CurrentAnalysisFlag(CURRENT_EMPTY_ANALYSIS)
{
  theAnalysisModel = new AnalysisModel();
}

BasicAnalysisBuilder::~BasicAnalysisBuilder()
{
    this->wipe();
    this->resetAll();
}

void BasicAnalysisBuilder::wipe()
{
    if (theStaticAnalysis != nullptr) {
        // theStaticAnalysis->clearAll();
        delete theStaticAnalysis;
        theStaticAnalysis = nullptr;
        resetStatic();
    }
    if (theTransientAnalysis != nullptr) {
        // theTransientAnalysis->clearAll();
        delete theTransientAnalysis;
        theTransientAnalysis = nullptr;
        resetTransient();
    }
    theVariableTimeStepTransientAnalysis = nullptr;
}

void BasicAnalysisBuilder::resetStatic()
{
    theAlgorithm = nullptr;
    theStaticIntegrator = nullptr;
    theSOE = nullptr;
    theNumberer = nullptr;
    theHandler = nullptr;
    theTest = nullptr;
//  theAnalysisModel = nullptr;
}

void BasicAnalysisBuilder::resetTransient()
{
    theAlgorithm = nullptr;
    theTransientIntegrator = nullptr;
    theSOE = nullptr;
    theNumberer = nullptr;
    theHandler = nullptr;
    theTest = nullptr;
//  theAnalysisModel = nullptr;
}

void BasicAnalysisBuilder::resetAll()
{
    theAlgorithm           = nullptr;
    theStaticIntegrator    = nullptr;
    theTransientIntegrator = nullptr;
    theSOE                 = nullptr;
    theNumberer            = nullptr;
    theHandler             = nullptr;
    theTest                = nullptr;

//  theAnalysisModel       = nullptr;
    theEigenSOE            = nullptr;
}

int
BasicAnalysisBuilder::domainChanged(void)
{
    int result = 0;

    Domain *domain = this->getDomain();
    int stamp = domain->hasDomainChanged();
    domainStamp = stamp;

    theAnalysisModel->clearAll();
    theHandler->clearAll();

    // now we invoke handle() on the constraint handler which
    // causes the creation of FE_Element and DOF_Group objects
    // and their addition to the AnalysisModel.
    result = theHandler->handle();
    if (result < 0) {
      opserr << "BasicAnalysisBuilder::domainChange() - ";
      opserr << "ConstraintHandler::handle() failed\n";
      return -1;
    }

    // we now invoke number() on the numberer which causes
    // equation numbers to be assigned to all the DOFs in the
    // AnalysisModel.
    result = theNumberer->numberDOF();
    if (result < 0) {
      opserr << "BasicAnalysisBuilder::domainChange() - ";
      opserr << "DOF_Numberer::numberDOF() failed\n";
      return -2;
    }

    result = theHandler->doneNumberingDOF();
    if (result < 0) {
      opserr << "BasicAnalysisBuilder::domainChange() - ";
      opserr << "ConstraintHandler::doneNumberingDOF() failed";
      return -2;
    }

    // we invoke setSize() on the LinearSOE which
    // causes that object to determine its size
    Graph &theGraph = theAnalysisModel->getDOFGraph();

    if (theSOE != nullptr) {
      result = theSOE->setSize(theGraph);
      if (result < 0) {
          opserr << "BasicAnalysisBuilder::domainChange() - ";
          opserr << "LinearSOE::setSize() failed";
          return -3;
      }
    }

    theAnalysisModel->clearDOFGraph();

    // finally we invoke domainChanged on the Integrator and Algorithm
    // objects .. informing them that the model has changed
    if (theStaticIntegrator != nullptr) {
      result = theStaticIntegrator->domainChanged();
      if (result < 0) {
        opserr << "BasicAnalysisBuilder::domainChange() - ";
        opserr << "Integrator::domainChanged() failed";
        return -4;
      }
    }
    if (theTransientIntegrator != nullptr) {
      result = theTransientIntegrator->domainChanged();
      if (result < 0) {
        opserr << "BasicAnalysisBuilder::domainChange() - ";
        opserr << "Integrator::domainChanged() failed";
        return -4;
      }
    }

    result = theAlgorithm->domainChanged();
    if (result < 0) {
      opserr << "StaticAnalysis::setAlgorithm() - ";
      opserr << "Algorithm::domainChanged() failed";
      return -5;
    }

    return 0;
}

int
BasicAnalysisBuilder::analyze(int num_steps, double size_steps)
{

  switch (this->CurrentAnalysisFlag) {

    case CURRENT_STATIC_ANALYSIS:
//    return this->analyzeStatic(num_steps);
      return theStaticAnalysis->analyze(num_steps);
      break;

    case CURRENT_TRANSIENT_ANALYSIS:
      // TODO: Set global timestep variable
      ops_Dt = size_steps;
//    return this->analyzeTransient(num_steps, size_steps);
      return theTransientAnalysis->analyze(num_steps, size_steps);
      break;

    default:
      opserr << G3_ERROR_PROMPT << "No Analysis type has been specified \n";
      return -1;
  }
}

int
BasicAnalysisBuilder::analyzeStatic(int numSteps)
{
  int result = 0;

  for (int i=0; i<numSteps; i++) {

      result = theAnalysisModel->analysisStep();

      if (result < 0) {
        opserr << "StaticAnalysis::analyze - the AnalysisModel failed";
        opserr << " at step: " << i << " with domain at load factor ";
        opserr << theDomain->getCurrentTime() << endln;
        theDomain->revertToLastCommit();
        return -2;
      }

      // Check for change in Domain since last step. As a change can
      // occur in a commit() in a domaindecomp with load balancing
      // this must now be inside the loop
      int stamp = theDomain->hasDomainChanged();

      if (stamp != domainStamp) {
        domainStamp = stamp;

        result = this->domainChanged();
        if (result < 0) {
          opserr << "BasicAnalysisBuilder::analyzeStatic - domainChanged failed";
          opserr << " at step " << i << " of " << numSteps << endln;
          return -1;
        }
      }

      result = theStaticIntegrator->newStep();
      if (result < 0) {
        opserr << "StaticAnalysis::analyze - the Integrator failed";
        opserr << " at step: " << i << " with domain at load factor ";
        opserr << theDomain->getCurrentTime() << endln;
        theDomain->revertToLastCommit();
        theStaticIntegrator->revertToLastStep();
        return -2;
      }

      result = theAlgorithm->solveCurrentStep();
      if (result < 0) {
        // opserr << "StaticAnalysis::analyze() - the Algorithm failed";
        // opserr << " at step: " << i << " with domain at load factor ";
        // opserr << theDomain->getCurrentTime() << endln;
        theDomain->revertToLastCommit();
        theStaticIntegrator->revertToLastStep();
        return -3;
      }

      result = theStaticIntegrator->commit();
      if (result < 0) {
        opserr << "StaticAnalysis::analyze - ";
        opserr << "the Integrator failed to commit";
        opserr << " at step: " << i << " with domain at load factor ";
        opserr << theDomain->getCurrentTime() << endln;

        theDomain->revertToLastCommit();
        theStaticIntegrator->revertToLastStep();
        return -4;
      }
  }

  return 0;
}

int
BasicAnalysisBuilder::analyzeTransient(int num_steps, double size_steps)
{
  return -1;
}

void BasicAnalysisBuilder::set(ConstraintHandler* obj) {

    if (obj == nullptr)
      return;

#if 0
    if (theHandler != nullptr) {
        // TODO: this needs to return a failure code.
        opserr << G3_WARN_PROMPT << "The handler can only be set once\n";
        return;
    }
#endif

    theHandler = obj;
}

void BasicAnalysisBuilder::set(DOF_Numberer* obj) {

    if (obj == nullptr)
      return;

    // if (theNumberer != nullptr) {
    //     opserr << "The numberer can only be set once for one analysis\n";
    //     return;
    // }

    theNumberer = obj;
    if (theStaticAnalysis != nullptr)
        theStaticAnalysis->setNumberer(*obj);

    if (theTransientAnalysis != nullptr)
        theTransientAnalysis->setNumberer(*obj);
}

void
BasicAnalysisBuilder::set(EquiSolnAlgo* obj)
{

    if (obj == nullptr)
      return;

    // if (theAlgorithm != nullptr) {
    //     opserr << "The algorithm can only be set once for one analysis\n";
    //     return;
    // }

    theAlgorithm = obj;
    if (theStaticAnalysis != nullptr)
        theStaticAnalysis->setAlgorithm(*obj);

    if (theTransientAnalysis != nullptr)
        theTransientAnalysis->setAlgorithm(*obj);
}

void
BasicAnalysisBuilder::set(LinearSOE* obj)
{
    if (obj == nullptr)
      return;
/*
    if (theSOE != nullptr)
      delete theSOE;
*/
    theSOE = obj;

    // NOTE: `setLinearSOE` will free any pre-existing LinearSOE.
    // maybe change this
    if (theStaticAnalysis != nullptr)
      theStaticAnalysis->setLinearSOE(*obj);

    if (theTransientAnalysis != nullptr)
      theTransientAnalysis->setLinearSOE(*obj);

#ifdef _PARALLEL_PROCESSING
    if (theStaticAnalysis != nullptr || theTransientAnalysis != nullptr) {
      SubdomainIter &theSubdomains = theDomain->getSubdomains();
      Subdomain *theSub;
      while ((theSub = theSubdomains()) != nullptr) {
	theSub->setAnalysisLinearSOE(*theSOE);
      }
    }
#endif

}

LinearSOE*
BasicAnalysisBuilder::getLinearSOE() {
  return theSOE;
}

// TODO: set(Integrator) is hideous
void
BasicAnalysisBuilder::set(Integrator* obj, int isstatic) {

    if (obj == nullptr)
      return;

    if (isstatic == 1) {
        // if (theStaticAnalysis != nullptr && theStaticIntegrator != nullptr)
        //       opserr << "WARNING - unexpected state.\n";

        theStaticIntegrator = dynamic_cast<StaticIntegrator*>(obj);
        if (theStaticIntegrator != nullptr) {
            if (theStaticAnalysis != nullptr) {
                theStaticAnalysis->setIntegrator(*theStaticIntegrator);
                return;
            }
        }
    }

    else {
        // if (theTransientIntegrator != nullptr)
        //     opserr << "WARNING - unexpected state\n";

        theTransientIntegrator = dynamic_cast<TransientIntegrator*>(obj);
        if (theTransientIntegrator != nullptr) {
            if (theTransientAnalysis != nullptr) {
                theTransientAnalysis->setIntegrator(*theTransientIntegrator);
                return;
            }
        }
    }
}


void
BasicAnalysisBuilder::set(ConvergenceTest* obj)
{

    if (obj == nullptr)
      return;

    theTest = obj;
    if (theStaticAnalysis != nullptr)
      theStaticAnalysis->setConvergenceTest(*obj);

    if (theTransientAnalysis != nullptr)
      theTransientAnalysis->setConvergenceTest(*obj);
}


void BasicAnalysisBuilder::newStaticAnalysis()
{
    // this->wipe();
    assert(theDomain != nullptr);

    if (theStaticAnalysis != nullptr) {
      delete theStaticAnalysis;
      theStaticAnalysis = nullptr;
    }

    if (theAnalysisModel == nullptr) {
      theAnalysisModel = new AnalysisModel();
    }

    if (theTest == nullptr) {
      theTest = new CTestNormUnbalance(1.0e-6,25,0);
    }

    if (theAlgorithm == nullptr) {
      theAlgorithm = new NewtonRaphson(*theTest);
    }

    if (theHandler == nullptr) {
      opserr << "WARNING analysis Static - no ConstraintHandler yet specified, \n";
      opserr << " PlainHandler default will be used\n";
      theHandler = new PlainHandler();
    }

    if (theNumberer == nullptr) {
      RCM *theRCM = new RCM(false);
      theNumberer = new DOF_Numberer(*theRCM);
    }

    if (theStaticIntegrator == nullptr) {
      // TODO
      opserr << "WARNING analysis Static - no Integrator specified, \n";
      //opserr << " StaticIntegrator default will be used\n";
      opserr << " LoadControl default will be used\n";
      theStaticIntegrator = new LoadControl(1, 1, 1, 1);
    }

    if (theSOE == nullptr) {
      // TODO: CHANGE TO MORE GENERAL SOE
        ProfileSPDLinSolver *theSolver;
        theSolver = new ProfileSPDLinDirectSolver();
        theSOE = new ProfileSPDLinSOE(*theSolver);
    }

    theStaticAnalysis = new StaticAnalysis(*theDomain,*theHandler,*theNumberer,*theAnalysisModel,
                                           *theAlgorithm,*theSOE,*theStaticIntegrator,theTest);

    // this->resetStatic();
}

int
BasicAnalysisBuilder::setStaticAnalysis()
{
  if (theStaticAnalysis == nullptr)
    this->newStaticAnalysis();

  this->CurrentAnalysisFlag = CURRENT_STATIC_ANALYSIS;

  return 0;
}

int
BasicAnalysisBuilder::setTransientAnalysis()
{
  if (theTransientAnalysis == nullptr)
    this->newTransientAnalysis();

  this->CurrentAnalysisFlag = CURRENT_TRANSIENT_ANALYSIS;

  return 1;
}

int
BasicAnalysisBuilder::newTransientAnalysis()
{
    // this->wipe();
    assert(theDomain != nullptr);

    if (theTransientAnalysis != nullptr) {
      delete theTransientAnalysis;
      theTransientAnalysis = nullptr;
    }

    if (theAnalysisModel == nullptr) {
        theAnalysisModel = new AnalysisModel();
    }

    if (theTest == nullptr) {
        theTest = new CTestNormUnbalance(1.0e-6,25,0);
    }

    if (theAlgorithm == nullptr) {
        theAlgorithm = new NewtonRaphson(*theTest);
    }

    if (theHandler == nullptr) {
        opserr << "WARNING analysis Transient dt tFinal - no ConstraintHandler\n";
        opserr << " yet specified, PlainHandler default will be used\n";
        theHandler = new PlainHandler();
    }

    if (theNumberer == nullptr) {
        RCM *theRCM = new RCM(false);
        theNumberer = new DOF_Numberer(*theRCM);
    }

    if (theTransientIntegrator == nullptr) {
        theTransientIntegrator = new Newmark(0.5,0.25);
    }

    if (theSOE == nullptr) {
#if 0
        opserr << "WARNING analysis Transient dt tFinal - no LinearSOE specified, \n";
        opserr << " ProfileSPDLinSOE default will be used\n";
#endif
        // TODO: Change to general default
        ProfileSPDLinSolver *theSolver;
        theSolver = new ProfileSPDLinDirectSolver();
        theSOE = new ProfileSPDLinSOE(*theSolver);
    }

    theTransientAnalysis=new DirectIntegrationAnalysis(*theDomain,*theHandler,*theNumberer,
                                                       *theAnalysisModel,*theAlgorithm,
                                                       *theSOE,*theTransientIntegrator,
                                                       theTest);
    // this->resetTransient();
    return 1;
}


void
BasicAnalysisBuilder::newEigenAnalysis(int typeSolver, double shift)
{

    assert(theAnalysisModel != nullptr);

    LoadControl theIntegrator(1, 1, 1, 1);

    if (theHandler == nullptr) {
      // TODO: Make temporary
      theHandler = new TransformationConstraintHandler();
    }

    if (theNumberer == nullptr) {
      // TODO: Make temporary
        theNumberer = new DOF_Numberer(*(new RCM(false)));
    }
    // TODO: FREE!!
    if (theSOE == nullptr) {
      // TODO: CHANGE TO MORE GENERAL SOE
        theSOE = new ProfileSPDLinSOE(*(new ProfileSPDLinDirectSolver()));
    }

    theAnalysisModel->setLinks(*theDomain, *theHandler);
    theHandler->setLinks(*theDomain, *theAnalysisModel, theIntegrator);
    theNumberer->setLinks(*theAnalysisModel);

    // create a new eigen system and solver
    if (theEigenSOE != nullptr) {
      if (theEigenSOE->getClassTag() != typeSolver) {
        // TODO
        //        delete theEigenSOE;
        theEigenSOE = nullptr;
      }
    }

    if (theEigenSOE == nullptr) {
      if (typeSolver == EigenSOE_TAGS_SymBandEigenSOE) {
        SymBandEigenSolver *theEigenSolver = new SymBandEigenSolver();
        theEigenSOE = new SymBandEigenSOE(*theEigenSolver, *theAnalysisModel);

      } else if (typeSolver == EigenSOE_TAGS_FullGenEigenSOE) {
          FullGenEigenSolver *theEigenSolver = new FullGenEigenSolver();
          theEigenSOE = new FullGenEigenSOE(*theEigenSolver, *theAnalysisModel);

      } else {
          theEigenSOE = new ArpackSOE(shift);
      }

      //
      // set the eigen soe in the system
      //
      theEigenSOE->setLinks(*theAnalysisModel);
      theEigenSOE->setLinearSOE(*theSOE);
    } // theEigenSOE == 0
}

int
BasicAnalysisBuilder::eigen(int numMode, bool generalized, bool findSmallest)
{
    // TODO: merge with newEigenAnalysis

    assert(theAnalysisModel != nullptr);
    assert(     theEigenSOE != nullptr);

    int result = 0;
    Domain *the_Domain = this->getDomain();

    // for parallel processing, want all analysis doing an eigenvalue analysis
    result = theAnalysisModel->eigenAnalysis(numMode, generalized, findSmallest);

    int stamp = the_Domain->hasDomainChanged();

    if (stamp != domainStamp) {
      domainStamp = stamp;
//    result = this->domainChanged();

      theAnalysisModel->clearAll();
      theHandler->clearAll();

      // Now invoke handle() on the constraint handler which
      // causes the creation of FE_Element and DOF_Group objects
      // and their addition to the AnalysisModel.
      result = theHandler->handle();
      if (result < 0) {
        opserr << "BasicAnalysisBuilder::eigen - ConstraintHandler::handle failed\n";
        return -1;
      }
      // Now invoke number() on the numberer which causes
      // equation numbers to be assigned to all the DOFs in the
      // AnalysisModel.
      result = theNumberer->numberDOF();
      if (result < 0) {
        opserr << "BasicAnalysisBuilder::eigen() - ";
        opserr << "DOF_Numberer::numberDOF() failed\n";
        return -2;
      }

      result = theHandler->doneNumberingDOF();
      if (result < 0) {
        opserr << "BasicAnalysisBuilder::eigen() - ";
        opserr << "ConstraintHandler::doneNumberingDOF() failed\n";
        return -2;
      }
      Graph &theGraph = theAnalysisModel->getDOFGraph();

      result = theSOE->setSize(theGraph);
      if (result < 0) {
          opserr << "BasicAnalysisBuilder::eigen() - ";
          opserr << "LinearSOE::setSize() failed\n";
          return -3;
      }

      result = theEigenSOE->setSize(theGraph);
      if (result < 0) {
	opserr << "BasicAnalysisBuilder::eigen() - ";
	opserr << "EigenSOE::setSize() failed\n";
	return -3;
      }

      theAnalysisModel->clearDOFGraph();

      if (result < 0) {
        opserr << "BasicAnalysisBuilder::eigen() - domainChanged failed\n";
	return -1;
      }	
    }

    //
    // zero A and M
    //
    theEigenSOE->zeroA();
    theEigenSOE->zeroM();

    //
    // form K
    //
    FE_EleIter &theEles = theAnalysisModel->getFEs();
    FE_Element *elePtr;
    while ((elePtr = theEles()) != nullptr) {
      elePtr->zeroTangent();
      elePtr->addKtToTang(1.0);
      if (theEigenSOE->addA(elePtr->getTangent(0), elePtr->getID()) < 0) {
	opserr << "WARNING DirectIntegrationAnalysis::eigen() -";
	opserr << " failed in addA for ID " << elePtr->getID();	
	result = -2;
      }
    }

    //
    // if generalized is true, form M
    //
    if (generalized == true) {
      FE_EleIter &theEles2 = theAnalysisModel->getFEs();
      while ((elePtr = theEles2()) != nullptr) {
	elePtr->zeroTangent();
	elePtr->addMtoTang(1.0);
	if (theEigenSOE->addM(elePtr->getTangent(0), elePtr->getID()) < 0) {
	  opserr << "WARNING BasicAnalysisBuilder::eigen() -";
	  opserr << " failed in addA for ID " << elePtr->getID() << "\n";
	  result = -2;
	}
      }

      DOF_Group *dofPtr;
      DOF_GrpIter &theDofs = theAnalysisModel->getDOFs();
      while ((dofPtr = theDofs()) != nullptr) {
	dofPtr->zeroTangent();
	dofPtr->addMtoTang(1.0);
	if (theEigenSOE->addM(dofPtr->getTangent(0), dofPtr->getID()) < 0) {
	  opserr << "WARNING BasicAnalysisBuilder::eigen() -";
	  opserr << " failed in addM for ID " << dofPtr->getID() << "\n";
	  result = -3;
	}
      }
    }

    //
    // solve for the eigen values & vectors
    //
    if (theEigenSOE->solve(numMode, generalized, findSmallest) < 0) {
	opserr << "WARNING BasicAnalysisBuilder::eigen() - EigenSOE failed in solve()\n";
	return -4;
    }
	
    //
    // now set the eigenvalues and eigenvectors in the model
    //
    theAnalysisModel->setNumEigenvectors(numMode);
    Vector theEigenvalues(numMode);
    for (int i = 1; i <= numMode; i++) {
      theEigenvalues[i-1] = theEigenSOE->getEigenvalue(i);
      theAnalysisModel->setEigenvector(i, theEigenSOE->getEigenvector(i));
    }
    theAnalysisModel->setEigenvalues(theEigenvalues);
    this->numEigen = numMode;

    return 0;
}

Domain*
BasicAnalysisBuilder::getDomain()
{
  return theDomain;
}

EquiSolnAlgo*
BasicAnalysisBuilder::getAlgorithm()
{
    return theAlgorithm;
}

StaticIntegrator*
BasicAnalysisBuilder::getStaticIntegrator() {
    if (theStaticAnalysis != nullptr) {
        return theStaticAnalysis->getIntegrator();
    }
    return nullptr;
}

TransientIntegrator*
BasicAnalysisBuilder::getTransientIntegrator() {
    if (theTransientAnalysis != nullptr) {
        return theTransientAnalysis->getIntegrator();
    }
    return nullptr;
}

ConvergenceTest*
BasicAnalysisBuilder::getConvergenceTest()
{
    if (theStaticAnalysis != nullptr) {
        return theStaticAnalysis->getConvergenceTest();

    } else if (theTransientAnalysis != nullptr) {
        return theTransientAnalysis->getConvergenceTest();
    }

    return theTest;
}

#if 0

// static int ops_printEle(OPS_Stream &output, PyObject *args);
// static int ops_printNode(OPS_Stream &output, PyObject *args);
// static int ops_printAlgo(OPS_Stream &output, PyObject *args);
// static int ops_printInteg(OPS_Stream &output, PyObject *args);

PyObject *ops_printModel(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    // file stream
    FileStream outputFile;
    OPS_Stream *output = &opserr;

    // domain
    Domain& theDomain = *(OPS_GetDomain());

    // just printModel()
    int numData = OPS_GetNumRemainingInputArgs();
    if (numData == 0) {
        opserr << theDomain;
        Py_INCREF(Py_None);
        return Py_None;
    }

    // printModel()
    int res = 0;
    while(numData > 0) {
        std::string type = OPS_GetString();
        if (type=="-ele" || type=="ele") {
            res = ops_printEle(*output,args);
            break;
        } else if (type=="-node" || type=="node") {
            res = ops_printNode(*output,args);
            break;
        } else if (type=="-integrator" || type=="integrator") {
            res = ops_printInteg(*output,args);
            break;
        } else if (type=="-algorithm" || type=="algorithm") {
            res = ops_printAlgo(*output,args);
            break;
        } else {
            // open file
            if (type=="-file" || type=="file") {
                type = OPS_GetString();
                numData = OPS_GetNumRemainingInputArgs();
            }
            if (outputFile.setFile(type.c_str(), openMode::APPEND) != 0) {
                PyErr_SetString(PyExc_RuntimeError,"failed to open file ");
                return NULL;
            }

            // just print(filename)
            if (numData == 0) {
                outputFile << theDomain;
                Py_INCREF(Py_None);
                return Py_None;
            }
            output = &outputFile;
        }
    }

    if (res < 0) {
        PyErr_SetString(PyExc_RuntimeError,"failed to print ");
        return NULL;
    }

    // close the output file
    outputFile.close();

    Py_INCREF(Py_None);
    return Py_None;
}

int ops_printEle(OPS_Stream &output, PyObject *args)
{
    int flag = 0;

    Domain& theDomain = *(OPS_GetDomain());

    // just print(<filename>, 'ele')
    if (OPS_GetNumRemainingInputArgs() == 0) {
        ElementIter &theElements = theDomain.getElements();
        Element *theElement;
        while ((theElement = theElements()) != 0) {
            theElement->Print(output);
        }
        return 0;
    }

    // if 'print <filename> Element flag int <int int ..>' get the flag
    std::string type = OPS_GetString();
    if (type=="flag" || type=="-flag") {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING printModel(<filename>, 'ele', <'flag', int> no int specified \n";
            return -1;
        }
        int numData = 1;
        if (OPS_GetIntInput(&numData,&flag) < 0) return -1;
    } else {
        int numData = OPS_GetNumRemainingInputArgs();
        int numArgs = PyTuple_Size(args);
        OPS_ResetCommandLine(numArgs, numArgs-numData-1, args);
    }

    // now print the Elements with the specified flag, 0 by default
    int numEle = OPS_GetNumRemainingInputArgs();
    if (numEle == 0) {
        ElementIter &theElements = theDomain.getElements();
        Element *theElement;
        while ((theElement = theElements()) != 0) {
            theElement->Print(output, flag);
        }
        return 0;
    } else {

        // otherwise print out the specified nodes i j k .. with flag
        ID theEle(numEle);
        if (OPS_GetIntInput(&numEle, &theEle(0)) < 0) return -1;
        theDomain.Print(output, 0, &theEle, flag);
    }

    return 0;
}

int ops_printNode(OPS_Stream &output, PyObject *args)
{
    int flag = 0;

    Domain& theDomain = *(OPS_GetDomain());

    // just print(<filename>, 'node')
    if (OPS_GetNumRemainingInputArgs() == 0) {
        NodeIter &theNodes = theDomain.getNodes();
        Node *theNode;
        while((theNode = theNodes()) != 0) {
            theNode->Print(output);
        }
        return 0;
    }

    // if 'print <filename> node flag int <int int ..>' get the flag
    std::string type = OPS_GetString();
    if (type=="flag" || type=="-flag") {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING printModel(<filename>, 'ele', <'flag', int> no int specified \n";
            return -1;
        }
        int numData = 1;
        if (OPS_GetIntInput(&numData,&flag) < 0) return -1;
    } else {
        int numData = OPS_GetNumRemainingInputArgs();
        int numArgs = PyTuple_Size(args);
        OPS_ResetCommandLine(numArgs, numArgs-numData-1, args);
    }

    // now print the nodes with the specified flag, 0 by default
    int numNode = OPS_GetNumRemainingInputArgs();
    if (numNode == 0) {
        NodeIter &theNodes = theDomain.getNodes();
        Node *theNode;
        while ((theNode = theNodes()) != 0) {
            theNode->Print(output, flag);
        }
        return 0;
    } else {

        // otherwise print out the specified nodes i j k .. with flag
        ID theNode(numNode);
        if (OPS_GetIntInput(&numNode, &theNode(0)) < 0) return -1;
        theDomain.Print(output, &theNode, 0, flag);
    }

    return 0;
}

int ops_printAlgo(OPS_Stream &output, PyObject *args)
{
    EquiSolnAlgo* theAlgorithm = anaBuilder.getAlgorithm();
    if (theAlgorithm == 0) return 0;

    // if just 'print <filename> algorithm'- no flag
    if (OPS_GetNumRemainingInputArgs() == 0) {
        theAlgorithm->Print(output);
        return 0;
    }

    // if 'print <filename> Algorithm flag' get the flag
    int flag;
    int numData = 1;
    if (OPS_GetIntInput(&numData, &flag) < 0) return -1;
    theAlgorithm->Print(output,flag);

    return 0;
}

int ops_printInteg(OPS_Stream &output, PyObject *args)
{
    StaticIntegrator* theStaticIntegrator=anaBuilder.getStaticIntegrator();
    TransientIntegrator* theTransientIntegrator=anaBuilder.getTransientIntegrator();

    if (theStaticIntegrator == 0 && theTransientIntegrator==0)
        return 0;

    IncrementalIntegrator *theIntegrator;
    if (theStaticIntegrator != 0)
        theIntegrator = theStaticIntegrator;
    else
        theIntegrator = theTransientIntegrator;

    // if just 'print <filename> algorithm'- no flag
    if (OPS_GetNumRemainingInputArgs() == 0) {
        theIntegrator->Print(output);
        return 0;
    }

    // if 'print <filename> Algorithm flag' get the flag
    int flag;
    int numData = 1;
    if (OPS_GetIntInput(&numData, &flag) < 0) return -1;
    theIntegrator->Print(output,flag);

    return 0;
}

PyObject *ops_wipeAnalysis(PyObject *self, PyObject *args)
{
    anaBuilder.wipe();
    anaBuilder.resetAll();

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *ops_wipeModel(PyObject *self, PyObject *args)
{
    anaBuilder.wipe();
    anaBuilder.resetAll();

    // wipe domain
    // Domain* theDomain = OPS_GetDomain();
    theDomain->clearAll();

    // wipe uniaxial material
    OPS_clearAllUniaxialMaterial();
    OPS_clearAllNDMaterial();

    // wipe sections
    OPS_clearAllSectionForceDeformation();
    OPS_clearAllSectionRepres();

    // wipe time series
    OPS_clearAllTimeSeries();

    // wipe GeomTransf
    OPS_ClearAllCrdTransf();

    // wipe BeamIntegration
    OPS_clearAllBeamIntegrationRule();

    ops_Dt = 0.0;

    Py_INCREF(Py_None);
    return Py_None;
}


PyObject *ops_nodeDisp(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    int numData = OPS_GetNumRemainingInputArgs();
    if (numData < 1) {
        PyErr_SetString(PyExc_RuntimeError,"WARNING want - nodeDisp nodeTag? <dof?>");
        return NULL;
    }

    // get inputs
    int data[2] = {0,0};
    if (OPS_GetIntInput(&numData,&data[0]) < 0) return NULL;
    data[1]--;

    // get nodal response
    Domain& theDomain = *(OPS_GetDomain());
    const Vector *nodalResponse = theDomain.getNodeResponse(data[0], Disp);
    if (nodalResponse == 0) {
        PyErr_SetString(PyExc_RuntimeError,"failed to get nodal response");
        return NULL;
    }

    // get disp
    int size = nodalResponse->Size();
    if (data[1] >= 0) {
        if (data[1] >= size) {
            PyErr_SetString(PyExc_RuntimeError,
                            "WARNING nodeDisp nodeTag? dof? - dofTag? too large\n");
            return NULL;
        }
        double value = (*nodalResponse)(data[1]);
        return Py_BuildValue("d", value);
    }

    // get list
    PyObject* theList = PyList_New(0);
    if (theList == 0) {
        PyErr_SetString(PyExc_RuntimeError,"failed to create disp list");
        return NULL;
    }

    for (int i=0; i<size; i++) {
        if (PyList_Append(theList,Py_BuildValue("d",(*nodalResponse)(i))) < 0) {
            PyErr_SetString(PyExc_RuntimeError,"failed to create disp list");
            return NULL;
        }
    }

    return theList;
}

PyObject *ops_getLoadFactor(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    int numData = OPS_GetNumRemainingInputArgs();
    if (numData < 1) {
        PyErr_SetString(PyExc_RuntimeError,"want patternTag");
        return NULL;
    }

    // get inputs
    int patternTag;
    numData = 1;
    if (OPS_GetIntInput(&numData,&patternTag) < 0) return NULL;

    // get load pattern
    Domain& theDomain = *(OPS_GetDomain());
    LoadPattern* thePattern = theDomain.getLoadPattern(patternTag);
    if (thePattern == 0) {
        PyErr_SetString(PyExc_RuntimeError,"the load pattern is not found");
        return NULL;
    }

    // get load factor
    double value = thePattern->getLoadFactor();
    return Py_BuildValue("d", value);
}

PyObject *ops_setLoadConst(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    int numData = OPS_GetNumRemainingInputArgs();
    if (numData < 2) {
        Py_INCREF(Py_None);
        return Py_None;
    }

    // Domain* theDomain = OPS_GetDomain();
    theDomain->setLoadConstant();

    std::string type = OPS_GetString();

    if (type == "-time") {
        double newTime;
        numData = 1;
        if (OPS_GetDoubleInput(&numData,&newTime) < 0) return NULL;
        theDomain->setCurrentTime(newTime);
        theDomain->setCommittedTime(newTime);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *ops_eigenAnalysis(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    // check inputs
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 1) {
        PyErr_SetString(PyExc_RuntimeError,"ERROR eigen('type',numModes)");
        return NULL;
    }

    // get parameters
    bool generalizedAlgo = true; // 0 - frequency/generalized (default),
                                 // 1 - standard,
                                 // 2 - buckling
    int typeSolver = EigenSOE_TAGS_ArpackSOE;
    double shift = 0.0;
    bool findSmallest = true;
    this->numEigen = 0;

    // check type of eigenvalue analysis
    if (numArgs >1) {
        std::string type = OPS_GetString();
        if (type=="frequency"||type=="-frenquency"||type=="generalized"||type=="-generalized") {
            generalizedAlgo = true;
        } else if (type=="standard"||type=="-standard") {
            generalizedAlgo = false;
        } else if (type=="-findLargest") {
            findSmallest = false;
        } else if (type=="genBandArpack"||type=="--genBandArpack"||
                  type=="genBandArpackEigen"||type=="-genBandArpackEigen") {
            typeSolver = EigenSOE_TAGS_ArpackSOE;
        } else if (type=="symmBandLapack"||type=="-symmBandLapack"||
                  type=="symmBandLapackEigen"||type=="-symmBandLapackEigen") {
            typeSolver = EigenSOE_TAGS_SymBandEigenSOE;
        } else if (type=="fullGenLapack"||type=="-fullGenLapack"||
                  type=="fullGenLapackEigen"||type=="-fullGenLapackEigen") {
            typeSolver = EigenSOE_TAGS_FullGenEigenSOE;
        } else {
            PyErr_SetString(PyExc_RuntimeError,"eigen - unknown option specified");
            return NULL;
        }
    }

    int numData = 1;
    if (OPS_GetIntInput(&numData,&numEigen) < 0) return NULL;

    anaBuilder.newEigenAnalysis(typeSolver,shift);

    // create a transient analysis if no analysis exists

    StaticAnalysis* theStaticAnalysis = anaBuilder.getStaticAnalysis();
    DirectIntegrationAnalysis* theTransientAnalysis = anaBuilder.getTransientAnalysis();
    if (theStaticAnalysis == 0 && theTransientAnalysis == 0) {
        anaBuilder.newTransientAnalysis();
    }

    // call analysis
    int res = 0;
    if (theStaticAnalysis != 0) {
        res = theStaticAnalysis->eigen(numEigen,generalizedAlgo,findSmallest);
    } else if (theTransientAnalysis != 0) {
        res = theTransientAnalysis->eigen(numEigen,generalizedAlgo,findSmallest);
    }

    // return
    PyObject* theList = PyList_New(0);
    if (res == 0) {
        const Vector& eigenvalues = theDomain->getEigenvalues();
        // get list
        if (theList == 0) {
            PyErr_SetString(PyExc_RuntimeError,"failed to create disp list");
            return NULL;
        }

        for (int i=0; i<numEigen; i++) {
            if (PyList_Append(theList,Py_BuildValue("d",eigenvalues(i))) < 0) {
                PyErr_SetString(PyExc_RuntimeError,"failed to create disp list");
                return NULL;
            }
        }
    }

    return theList;

}
#endif
