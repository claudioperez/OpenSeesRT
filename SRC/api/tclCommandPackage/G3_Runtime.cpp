// #include <stdio.h>

#include <unordered_map>
#include <string>
#include <vector>
#include "G3_Runtime.h"
// #include "G3Parse.h"

#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <LinearSOE.h>
#include <EigenSOE.h>
#include <DOF_Numberer.h>
#include <ConstraintHandler.h>
#include <StaticIntegrator.h>
#include <TransientIntegrator.h>
#include <EquiSolnAlgo.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>

// DEFAULTS
#include <AnalysisModel.h>
#include <ProfileSPDLinSolver.h>
#include <ProfileSPDLinDirectSolver.h>
#include <NewtonRaphson.h>
#include <TransformationConstraintHandler.h>
#include <CTestNormUnbalance.h>
#include <ProfileSPDLinSOE.h>
#include <PlainHandler.h>
#include <Newmark.h>
#include <RCM.h>
#include <LoadControl.h>

class StaticIntegrator;

#define G3Config_keyExists(conf, key) ((conf).find((key)) != (conf).end())

// typedef void*(G3_Parse)(Tcl_Interp*, int, const char **);
template <typename T>
using G3_Parse = T* (*)(G3_Runtime*, int, const char **);


// Wrap a function with signature
//  T* G3Parse_newT(G3_Runtime*, int argc, G3_Char** argv)
// so that it works with a std::vector<std::string>
// written: cmp
template<typename T, T* (*fn)(G3_Runtime*, int, G3_Char **)>
T* G3Object_newParsed(G3_Runtime *rt, std::vector<std::string> args) {
    std::vector<G3_Char *> cstrs;
    cstrs.reserve(args.size());
    for (auto &s : args) cstrs.push_back(const_cast<char *>(s.c_str()));
    return (*fn)(rt, cstrs.size()+1, cstrs.data()-1);
}

G3_Parse<ConvergenceTest>     G3Parse_newConvergenceTest;
// G3_Parse<TransientIntegrator> G3Parse_newTransientIntegrator;
// G3_Parse<StaticIntegrator>    G3Parse_newStaticIntegrator;
EquiSolnAlgo* G3Parse_newEquiSolnAlgo(G3_Runtime*, int, const char **);
TransientIntegrator* G3Parse_newTransientIntegrator(G3_Runtime*, int, const char**);
StaticIntegrator* G3Parse_newStaticIntegrator(G3_Runtime*, int, const char**);
LinearSOE* G3Parse_newLinearSOE(G3_Runtime*, int, const char**);



// template<typename T, std::string conf_key>



void *
G3_Runtime::newStaticAnalysis(G3_Config conf)
{
  StaticIntegrator* sintegrator = nullptr;

  // INTEGRATOR
  if (G3Config_keyExists(conf, "integrator"))
    sintegrator = 
      G3Object_newParsed<StaticIntegrator, G3Parse_newStaticIntegrator>(this, conf["integrator"]);
  else
    sintegrator = new LoadControl(1, 1, 1, 1);

  // CONVERGENCE TEST
  ConvergenceTest *test = new CTestNormUnbalance(1.0e-6,25,0);

  // ALGORITHM
  EquiSolnAlgo* the_algorithm;
  if (G3Config_keyExists(conf, "algorithm"))
    the_algorithm = 
      G3Object_newParsed<EquiSolnAlgo, G3Parse_newEquiSolnAlgo>(this, conf["algorithm"]);
  else
    the_algorithm = this->m_global_strategy.m_algorithm;
  if (the_algorithm == nullptr)
    the_algorithm = new NewtonRaphson(*test);
  else
    the_algorithm->setConvergenceTest(test);


  // NUMBERER
  DOF_Numberer* the_numberer = nullptr;
  /*
  if (G3Config_keyExists(conf, "numberer"))
    the_numberer = 
      G3Object_newParsed<DOF_Numeberer, G3Parse_newNumberer>(this, conf["numberer"]);
  else
    the_numberer = this->m_global_strategy.m_numberer;
  */
  if (the_numberer == nullptr)  {
    RCM *theRCM  = new RCM(false);
    the_numberer = new DOF_Numberer(*theRCM);
  }
  
  
  // CONSTRAINT HANDLER
  ConstraintHandler *the_handler = new TransformationConstraintHandler();


  // LINEAR SYSTEM
  LinearSOE* the_soe = nullptr;
  if (G3Config_keyExists(conf, "system"))
    the_soe = 
      G3Object_newParsed<LinearSOE, G3Parse_newLinearSOE>(this, conf["system"]);
  else
    the_soe = this->m_global_strategy.m_linear_soe;

  if (the_soe == nullptr) 
      the_soe = new ProfileSPDLinSOE(*new ProfileSPDLinDirectSolver());

  
  if (m_analysis_model == nullptr)
    m_analysis_model = new AnalysisModel();

  return  new StaticAnalysis(*m_domain,
                             *the_handler,
                             *the_numberer,
                             *m_analysis_model,
                             *the_algorithm,
                             *the_soe,
                             *sintegrator,
                             test);
}

void *
G3_Runtime::newTransientAnalysis(G3_Config conf)
{
  // NUMBERER
  DOF_Numberer* the_numberer = nullptr;
  /*
  if (G3Config_keyExists(conf, "numberer"))
    the_numberer = 
      G3Object_newParsed<DOF_Numeberer, G3Parse_newNumberer>(this, conf["numberer"]);
  else
    the_numberer = this->m_global_strategy.m_numberer;
  */
  if (the_numberer == nullptr)  {
    RCM *theRCM  = new RCM(false);
    the_numberer = new DOF_Numberer(*theRCM);
  }

  // CONSTRAINT HANDLER
  ConstraintHandler *the_handler = new TransformationConstraintHandler();

  // CONVERGENCE TEST
  ConvergenceTest *test = new CTestNormUnbalance(1.0e-6,25,0);

  // ALGORITHM
  EquiSolnAlgo* the_algorithm;
  if (G3Config_keyExists(conf, "algorithm"))
    the_algorithm = 
      G3Object_newParsed<EquiSolnAlgo, G3Parse_newEquiSolnAlgo>(this, conf["algorithm"]);
  else
    the_algorithm = this->m_global_strategy.m_algorithm;
  if (the_algorithm == nullptr)
    the_algorithm = new NewtonRaphson(*test);
  else
    the_algorithm->setConvergenceTest(test);


  // LINEAR SYSTEM
  LinearSOE* the_soe = nullptr;
  if (G3Config_keyExists(conf, "system"))
    the_soe = 
      G3Object_newParsed<LinearSOE, G3Parse_newLinearSOE>(this, conf["system"]);
  else
    the_soe = this->m_global_strategy.m_linear_soe;

  if (the_soe == nullptr) 
      the_soe = new ProfileSPDLinSOE(*new ProfileSPDLinDirectSolver());


  // ANALYSIS MODEL
  if (m_analysis_model == nullptr)
    m_analysis_model = new AnalysisModel();
 

  TransientIntegrator* tintegrator;
  if (G3Config_keyExists(conf, "integrator"))
    tintegrator = 
      G3Object_newParsed<TransientIntegrator, G3Parse_newTransientIntegrator>(this, conf["integrator"]);
  else
      tintegrator = new Newmark(0.5, 0.25);



  if (G3Config_keyExists(conf, "analysis")) {
    if (!conf["analysis"].empty() && (conf["analysis"][0] == "Variable"))
      return new VariableTimeStepDirectIntegrationAnalysis(
                                       *m_domain,
                                       *the_handler,
                                       *the_numberer,
                                       *m_analysis_model,
                                       *the_algorithm,
                                       *the_soe,
                                       *tintegrator,
                                       test);
  }
    return new DirectIntegrationAnalysis(
                                       *m_domain,
                                       *the_handler,
                                       *the_numberer,
                                       *m_analysis_model,
                                       *the_algorithm,
                                       *the_soe,
                                       *tintegrator,
                                       test);
}

