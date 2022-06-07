#include <tcl.h>
#include <g3_api.h>
#include <stdio.h>

#include <unordered_map>
#include <string>
#include <vector>
typedef std::unordered_map<std::string, std::vector<std::string>> G3_Config;

class Domain;
class ModelBuilder;

class AnalysisModel;
class ConstraintHandler;
class LinearSOE;
class EigenSOE;
class DOF_Numberer;
class ConvergenceTest;

class G3_Interpreter;
#define G3_Builder G3_Runtime

class G3_Runtime {
public:
  // newStaticAnalysis()

  Tcl_Interp     *m_interp;
// MODEL BUILDING
  TclBuilder     *m_builder = nullptr;

  Domain         *m_domain  = nullptr;

  bool            model_is_built=false;

// ANALYSIS
  AnalysisModel  *m_analysis_model     = nullptr;
  AnalysisModel **m_analysis_model_ptr = &m_analysis_model;

  LinearSOE      *m_sys_of_eqn = nullptr;


  /*
  G3_Config m_global_strategy = {
    {"numberer",   {"RCM"}},
    {"system",     {"ProfileSPD"}},
    {"test",       {""}},
  };
  */

  struct G3_Strategy {
    EquiSolnAlgo      *m_algorithm  = nullptr;
    ConvergenceTest   *m_convergence_test = nullptr;
    LinearSOE         *m_linear_sys = nullptr;
    EigenSOE          *m_eigen_sys  = nullptr;
    ConstraintHandler *m_handler    = nullptr;
    DOF_Numberer      *m_numberer   = nullptr;
    LinearSOE         *m_linear_soe = nullptr;
  } m_global_strategy;


  void *newStaticAnalysis(G3_Config);
  void *newTransientAnalysis(G3_Config);


// IO
  FILE* streams[3] = {stdin,stdout,stderr};
};




class G3_ParallelRuntime : public G3_Runtime {
  bool is_partitioned=false;
  int num_subdomains = 0;
  bool flag_MPID_SOE = false;
};
