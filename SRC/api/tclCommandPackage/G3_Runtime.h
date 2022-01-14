#include <g3_api.h>

class ModelBuilder;
class G3_Interpreter;

class G3_Runtime {
public:
  Tcl_Interp *interp;
  bool model_is_built=false;
  // G3_Interpreter  *getInterpreter(void);
  // G3_ModelBuilder *getModelBuilder(void);

private:
// G3_Interpreter *interp;
// ModelBuilder *builder;

};
