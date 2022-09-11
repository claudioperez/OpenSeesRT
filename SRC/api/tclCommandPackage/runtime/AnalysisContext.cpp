#include <vector>
#include <string>
#include <unordered_map>

#include <runtime/AnalysisContext.hpp>



int AnalysisContext::update(char* key, int argc, char **argv)
{
}


SuccessFlag
update(std::unordered_map<std::string, std::vector<std::string>> conf)
{
}

SuccessFlag
AnalysisContext::setStaticIntegrator(StaticIntegrator* x)
{
  m_static_integrator = x;
  return G3_OK;
}

StaticIntegrator*
AnalysisContext::getStaticIntegrator()
{
   return m_static_integrator;
}


SuccessFlag
AnalysisContext::setConvergenceTest(ConvergenceTest* x)
{
  m_convergence_test = x;
  return G3_OK;
}

ConvergenceTest*
AnalysisContext::getConvergenceTest()
{
   return m_convergence_test;
}


SuccessFlag
AnalysisContext::setEquiSolnAlgo(EquiSolnAlgo* x)
{
  m_algorithm = x;
  return G3_OK;
}

EquiSolnAlgo*
AnalysisContext::getEquiSolnAlgo()
{
   return m_algorithm;
}

