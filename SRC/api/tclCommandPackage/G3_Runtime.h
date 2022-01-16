#include <g3_api.h>

class Domain;
class ModelBuilder;
class G3_Interpreter;

class G3_Runtime {
public:
  Tcl_Interp   *m_interp;
  TclBuilder   *m_builder = nullptr;
  Domain       *m_domain  = nullptr;
  bool model_is_built=false;
  // G3_Interpreter  *getInterpreter(void);
  // G3_ModelBuilder *getModelBuilder(void);

private:
// G3_Interpreter *interp;
// ModelBuilder *builder;

};
