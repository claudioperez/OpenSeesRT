#include <string>
#include <unordered_map>

#include <InputAPI.h>

class ConvergenceTest;

ConvergenceTest*
TclDispatch_newConvergenceTest(G3_Runtime* rt, int argc, G3_Char** argv);

std::unordered_map<std::string, ConvergenceTest* (*)(G3_Runtime*, int, G3_Char**)> 
ConvergenceTestLibrary = {
  {"LoadControl",                TclDispatch_newConvergenceTest},
  {"NormDispAndUnbalance",       TclDispatch_newConvergenceTest},
  {"NormDispOrUnbalance",        TclDispatch_newConvergenceTest},
  {"FixedNumIter",               TclDispatch_newConvergenceTest},
  {"FixedNumIter",               TclDispatch_newConvergenceTest},
  {"NormUnbalance",              TclDispatch_newConvergenceTest},
  {"NormDispIncr",               TclDispatch_newConvergenceTest},
  {"NormDispAndUnbalance",       TclDispatch_newConvergenceTest},
  {"NormDispOrUnbalance",        TclDispatch_newConvergenceTest},
  {"EnergyIncr",                 TclDispatch_newConvergenceTest},
  {"RelativeNormUnbalance",      TclDispatch_newConvergenceTest},
  {"RelativeNormDispIncr",       TclDispatch_newConvergenceTest},
  {"RelativeEnergyIncr",         TclDispatch_newConvergenceTest},
  {"RelativeTotalNormDispIncr",  TclDispatch_newConvergenceTest},
};
