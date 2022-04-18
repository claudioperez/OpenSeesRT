#include <tcl.h>
#include <g3_api.h>
#include <stdio.h>

class Domain;
class ModelBuilder;
class AnalysisModel;
class G3_Interpreter;

class G3_Runtime {
public:
  Tcl_Interp     *m_interp;
// MODEL BUILDING
  TclBuilder     *m_builder = nullptr;
  Domain         *m_domain  = nullptr;
  bool model_is_built=false;
// ANALYSIS
  AnalysisModel  *m_analysis_model     = nullptr;
  AnalysisModel **m_analysis_model_ptr = &m_analysis_model;

  LinearSOE      *m_sys_of_eqn = nullptr;

// IO
  FILE* stdout = stdout;
  FILE* stderr = stderr;
};

class G3_ParallelRuntime : public G3_Runtime {
  bool is_partitioned=false;
  int num_subdomains = 0;
  bool flag_MPID_SOE = false;
};
