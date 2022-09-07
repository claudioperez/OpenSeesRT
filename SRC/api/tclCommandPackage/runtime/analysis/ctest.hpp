#include <string>
#include <unordered_map>

#include <InputAPI.h>

class ConvergenceTest;

ConvergenceTest*
RT_newConvergenceTest(G3_Runtime* rt, int argc, G3_Char** argv);

std::unordered_map<std::string, ConvergenceTest* (*)(G3_Runtime*, int, G3_Char**)> 
ConvergenceTestLibrary = {
  {"LoadControl",                RT_newConvergenceTest},
  {"NormDispAndUnbalance",       RT_newConvergenceTest},
  {"NormDispOrUnbalance",        RT_newConvergenceTest},
  {"FixedNumIter",               RT_newConvergenceTest},
  {"FixedNumIter",               RT_newConvergenceTest},
  {"NormUnbalance",              RT_newConvergenceTest},
  {"NormDispIncr",               RT_newConvergenceTest},
  {"NormDispAndUnbalance",       RT_newConvergenceTest},
  {"NormDispOrUnbalance",        RT_newConvergenceTest},
  {"EnergyIncr",                 RT_newConvergenceTest},
  {"RelativeNormUnbalance",      RT_newConvergenceTest},
  {"RelativeNormDispIncr",       RT_newConvergenceTest},
  {"RelativeEnergyIncr",         RT_newConvergenceTest},
  {"RelativeTotalNormDispIncr",  RT_newConvergenceTest},
};
