#include <vector>
#include <string>
#include <unordered_map>

#include <InputAPI.h>
#include "analysis/numberer.hpp"
#include "analysis/system.hpp"
#include "analysis/integrator.hpp"
#include "analysis/ctest.hpp"
#include "analysis/algorithm.hpp"
// #include "analysis/eigen_soe.hpp"


template <typename T, std::unordered_map<std::string, T* (*)(G3_Runtime*, int, G3_Char**)>*>
class ConfigurationObject {
  int set(G3_Builder *mb, int argc, G3_Char **argv) {
    return 1;
  }
};


class AnalysisContext {
  int update(char*, int, char **);
  int update(std::unordered_map<std::string, std::vector<std::string>> conf);

  private:
    DOF_Numberer* m_numberer;
    ConfigurationObject<DOF_Numberer, &DOF_NumbererLibrary> m_numberer_library;

    int numEigen = 0;
    // EigenSOE* m_eigen_soe;
    // ConfigurationObject<EigenSOE, &EigenSOE_Library> m_linear_soe_library;

    LinearSOE* m_linear_soe;
    ConfigurationObject<LinearSOE, &LinearSOE_Library> m_linear_soe_library;

    TransientIntegrator* m_transient_integrator;
    ConfigurationObject<TransientIntegrator, &TransientIntegratorLibrary> m_transient_integrator_library;

    StaticIntegrator* m_static_integrator;
    ConfigurationObject<StaticIntegrator, &StaticIntegratorLibrary> m_static_integrator_library;

    ConvergenceTest* m_convergence_test;
    ConfigurationObject<ConvergenceTest, &ConvergenceTestLibrary> m_convergence_test_library;

    EquiSolnAlgo* m_algorithm;
    ConfigurationObject<EquiSolnAlgo, &EquiSolnAlgoLibrary> m_algorithm_library;
};

