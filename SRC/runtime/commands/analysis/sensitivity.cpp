//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
#include <tcl.h>
#include <Logging.h>
#include <Parsing.h>
#include <Integrator.h>
#include <StaticIntegrator.h>
#include <TransientIntegrator.h>
#include <BasicAnalysisBuilder.h>

#if 1
int 
computeGradients(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
    BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);
    Integrator* theIntegrator = nullptr;

    if(builder->getStaticIntegrator() != nullptr) {
        theIntegrator = builder->getStaticIntegrator();

    } else if(builder->getTransientIntegrator() != nullptr) {
        theIntegrator = builder->getTransientIntegrator();
    }

    if (theIntegrator == nullptr) {
        opserr << "WARNING: No integrator is created\n";
        return TCL_ERROR;
    }

    if (theIntegrator->computeSensitivities() < 0) {
      opserr << "WARNING: failed to compute sensitivities\n";
      return TCL_ERROR;
    }
    
    return TCL_OK;
}
#endif


int
TclCommand_sensitivityAlgorithm(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
    BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);


    Integrator* theIntegrator = nullptr;

    if (builder->getStaticIntegrator() != nullptr) {
        theIntegrator = builder->getStaticIntegrator();

    } else if(builder->getTransientIntegrator() != nullptr) {
        theIntegrator = builder->getTransientIntegrator();
    }


    if (argc < 2) {
        opserr << "ERROR: Wrong number of parameters to sensitivity algorithm." << "\n";
        return TCL_ERROR;
    }

    if (theIntegrator == nullptr) {
        opserr << "The integrator needs to be instantiated before " << "\n"
               << " setting  sensitivity algorithm." << "\n";
        return -1;
    }


    // 1: compute at each step (default); 
    // 2: compute by command; 

    int analysisTypeTag = 1;
    if (strcmp(argv[1],"-computeAtEachStep") == 0)
        analysisTypeTag = 1;

    else if (strcmp(argv[1],"-computeByCommand") == 0)
        analysisTypeTag = 2;

    else {
        opserr << "Unknown sensitivity algorithm option: " << argv[1] << "\n";
        return TCL_ERROR;
    }

    theIntegrator->setComputeType(analysisTypeTag);
    theIntegrator->activateSensitivityKey();
    
    return TCL_OK;
}

