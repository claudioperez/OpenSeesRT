/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Written: Minjie Zhu, Claudio Perez
//
#include "BasicAnalysisBuilder.h"
#include <Domain.h>
#include <string>
#include <assert.h>

#include <G3_Logging.h>

#include <Element.h>
#include <ElementIter.h>
#include <Node.h>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>

#include <EquiSolnAlgo.h>
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
:theHandler(nullptr),theNumberer(nullptr),theAnalysisModel(nullptr),theAlgorithm(nullptr),
 theSOE(nullptr),theEigenSOE(nullptr),theStaticIntegrator(nullptr),theTransientIntegrator(nullptr),
 theTest(nullptr),theStaticAnalysis(nullptr),theTransientAnalysis(nullptr),
 theVariableTimeStepTransientAnalysis(nullptr)
{

}

BasicAnalysisBuilder::BasicAnalysisBuilder(Domain* domain)
:theHandler(nullptr),theNumberer(nullptr),theAnalysisModel(nullptr),theAlgorithm(nullptr),
 theSOE(nullptr),theEigenSOE(nullptr),theStaticIntegrator(nullptr),theTransientIntegrator(nullptr),
 theTest(nullptr),theStaticAnalysis(nullptr),theTransientAnalysis(nullptr),
 theVariableTimeStepTransientAnalysis(nullptr), theDomain(domain)
{

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
    theAnalysisModel = nullptr;
}

void BasicAnalysisBuilder::resetTransient()
{
    theAlgorithm = nullptr;
    theTransientIntegrator = nullptr;
    theSOE = nullptr;
    theNumberer = nullptr;
    theHandler = nullptr;
    theTest = nullptr;
    theAnalysisModel = nullptr;
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
    theAnalysisModel       = nullptr;
    theEigenSOE            = nullptr;
}

#include <stdio.h>
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

void BasicAnalysisBuilder::set(EquiSolnAlgo* obj) {

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

void BasicAnalysisBuilder::set(LinearSOE* obj) {
    if (obj == nullptr)
      return;

    theSOE = obj;
    // NOTE: `setLinearSOE` will free any pre-existing LinearSOE.
    // maybe change this
    if (theStaticAnalysis != nullptr)
      theStaticAnalysis->setLinearSOE(*obj);

    if (theTransientAnalysis != nullptr)
      theTransientAnalysis->setLinearSOE(*obj);

#ifdef _PARALLEL_PROCESSING
    if (theStaticAnalysis != 0 || theTransientAnalysis != 0) {
      SubdomainIter &theSubdomains = theDomain.getSubdomains();
      Subdomain *theSub;
      while ((theSub = theSubdomains()) != 0) {
	theSub->setAnalysisLinearSOE(*theSOE);
      }
    }
#endif

}

LinearSOE*
BasicAnalysisBuilder::getLinearSOE(int flag) {
  return theSOE;
}

void BasicAnalysisBuilder::set(Integrator* obj, int isstatic) {

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

    if (obj == 0) return;

    // if (theTest != 0) {
    //     opserr << "The test can only be set once for one analysis\n";
    //     return;
    // } else {
    //
    // }
    theTest = obj;
    if (theStaticAnalysis != 0) 
      theStaticAnalysis->setConvergenceTest(*obj);

    if (theTransientAnalysis != 0) 
      theTransientAnalysis->setConvergenceTest(*obj);
}

int
BasicAnalysisBuilder::setStaticAnalysis()
{
  if (theStaticAnalysis == nullptr)
    this->newStaticAnalysis();

  this->CurrentAnalysisFlag = CURRENT_STATIC_ANALYSIS;

  return 0;
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
        // opserr << "WARNING analysis Static - no Numberer specified, \n";
        // opserr << " RCM default will be used\n";
        RCM *theRCM = new RCM(false);        
        theNumberer = new DOF_Numberer(*theRCM);            
    }
    if (theStaticIntegrator == nullptr) {
        opserr << "WARNING analysis Static - no Integrator specified, \n";
        opserr << " StaticIntegrator default will be used\n";
        theStaticIntegrator = new LoadControl(1, 1, 1, 1);       
    }

    if (theSOE == nullptr) {
      // TODO: CHANGE TO MORE GENERAL SOE
#if 0
        opserr << "WARNING analysis Static - no LinearSOE specified, \n";
        opserr << " ProfileSPDLinSOE default will be used\n";
#endif
        ProfileSPDLinSolver *theSolver;
        theSolver = new ProfileSPDLinDirectSolver();
        theSOE = new ProfileSPDLinSOE(*theSolver);      
    }

    // Domain* theDomain = OPS_GetDomain();
    theStaticAnalysis = new StaticAnalysis(*theDomain,*theHandler,*theNumberer,*theAnalysisModel,
                                           *theAlgorithm,*theSOE,*theStaticIntegrator,theTest);

    // set eigen SOE
    if (theEigenSOE != 0) {
        theStaticAnalysis->setEigenSOE(*theEigenSOE);
    }

    // this->resetStatic();
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
        opserr << "WARNING analysis Transient - no Algorithm yet specified, \n";
        opserr << " NewtonRaphson default will be used\n";            

        theAlgorithm = new NewtonRaphson(*theTest); 
    }
    if (theHandler == nullptr) {
        opserr << "WARNING analysis Transient dt tFinal - no ConstraintHandler\n";
        opserr << " yet specified, PlainHandler default will be used\n";
        theHandler = new PlainHandler();       
    }
    if (theNumberer == nullptr) {
        opserr << "WARNING analysis Transient dt tFinal - no Numberer specified, \n";
        opserr << " RCM default will be used\n";
        RCM *theRCM = new RCM(false);        
        theNumberer = new DOF_Numberer(*theRCM);            
    }
    if (theTransientIntegrator == nullptr) {
        opserr << "WARNING analysis Transient dt tFinal - no Integrator specified, \n";
        opserr << " Newmark(.5,.25) default will be used\n";
        theTransientIntegrator = new Newmark(0.5,0.25);       
    }
    if (theSOE == nullptr) {
        opserr << "WARNING analysis Transient dt tFinal - no LinearSOE specified, \n";
        opserr << " ProfileSPDLinSOE default will be used\n";
        ProfileSPDLinSolver *theSolver;
        theSolver = new ProfileSPDLinDirectSolver();         
        theSOE = new ProfileSPDLinSOE(*theSolver);      
    }

    theTransientAnalysis=new DirectIntegrationAnalysis(*theDomain,*theHandler,*theNumberer,
                                                       *theAnalysisModel,*theAlgorithm,
                                                       *theSOE,*theTransientIntegrator,
                                                       theTest);

    // set eigen SOE
    if (theEigenSOE != nullptr) {
        if (theTransientAnalysis != 0) {
          theTransientAnalysis->setEigenSOE(*theEigenSOE);
        }
    }

    // this->resetTransient();

    return 1;
}


void BasicAnalysisBuilder::newEigenAnalysis(int typeSolver, double shift)
{

    if (theHandler == nullptr) {
      theHandler = new TransformationConstraintHandler();
    }

    // create a new eigen system and solver
    if (theEigenSOE != nullptr) {
        if (theEigenSOE->getClassTag() != typeSolver) {
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
        if(theStaticAnalysis == nullptr && theTransientAnalysis == nullptr) {
          // TODO: these are only created to initialize a StaticAnalysis
          // object, but is not required
          this->set(new CTestNormUnbalance(1.0e-6,25,0));
          this->set(new LoadControl(1, 1, 1, 1), true);

          this->newStaticAnalysis();
          // theStaticAnalysis->setEigenSOE(*theEigenSOE);
        }

        if (theStaticAnalysis != nullptr) {
            theStaticAnalysis->setEigenSOE(*theEigenSOE);

        } /* else if (theTransientAnalysis == nullptr){
          this->newTransientAnalysis();
        }*/

        if (theTransientAnalysis != nullptr) {
            theTransientAnalysis->setEigenSOE(*theEigenSOE);
        }

    } // theEigenSOE != 0
}

Domain* BasicAnalysisBuilder::getDomain()
{
  return theDomain;
}

EquiSolnAlgo* BasicAnalysisBuilder::getAlgorithm()
{
    if (theStaticAnalysis != 0) {
        return theStaticAnalysis->getAlgorithm();
    } else if (theTransientAnalysis != 0) {
        return theTransientAnalysis->getAlgorithm();
    }

    return 0;
}

StaticIntegrator* BasicAnalysisBuilder::getStaticIntegrator() {
    if (theStaticAnalysis != 0) {
        return theStaticAnalysis->getIntegrator();
    }
    return 0;
}

TransientIntegrator* BasicAnalysisBuilder::getTransientIntegrator() {
    if (theTransientAnalysis != 0) {
        return theTransientAnalysis->getIntegrator();
    }
    return 0;
}

ConvergenceTest* BasicAnalysisBuilder::getConvergenceTest()
{
    if (theStaticAnalysis != 0) {
        return theStaticAnalysis->getConvergenceTest();

    } else if (theTransientAnalysis != 0) {
        return theTransientAnalysis->getConvergenceTest();
    }

    return theTest;
}

#if 0 
extern LinearSOE* OPS_ParseSOECommand(const char *type);
extern DOF_Numberer* OPS_ParseNumbererCommand(const char *type);
extern EquiSolnAlgo* OPS_ParseAlgorithmCommand(const char *type);
extern Integrator* OPS_ParseIntegratorCommand(const char *type, int& isstatic);

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
            if (outputFile.setFile(type.c_str(), APPEND) != 0) {
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
    // Domain* theDomain = OPS_GetDomain();

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
