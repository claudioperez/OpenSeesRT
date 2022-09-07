#include <string>
#include <unordered_map>

#include <InputAPI.h>
#include <analysis/numberer/DOF_Numberer.h>
#include <analysis/numberer/PlainNumberer.h>
#include <graph/numberer/RCM.h>
#include <graph/numberer/AMDNumberer.h>

// #include <ParallelNumberer.h>

class G3_Runtime;
typedef DOF_Numberer* (Function)(G3_Runtime*, int, G3_Char**);
std::unordered_map<std::string, Function*> DOF_NumbererLibrary = {
  {"Plain",  (Function*)[](G3_Runtime*, int, G3_Char**)->DOF_Numberer*{return new PlainNumberer();}        },
  {"AMD",    (Function*)[](G3_Runtime*, int, G3_Char**)->DOF_Numberer*{return new DOF_Numberer(*(new AMD()));}},
  {"RCM",    (Function*)[](G3_Runtime*, int, G3_Char**)->DOF_Numberer*{return new DOF_Numberer(*(new RCM(false)));}}
};

