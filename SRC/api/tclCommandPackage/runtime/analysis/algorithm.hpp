#include <string>
#include <unordered_map>

#include <InputAPI.h>

class EquiSolnAlgo;
EquiSolnAlgo* G3Parse_newEquiSolnAlgo(G3_Runtime*, int, G3_Char**);
EquiSolnAlgo* G3Parse_newLinearAlgorithm(G3_Runtime*, int, G3_Char**);
EquiSolnAlgo* G3Parse_newSecantNewtonAlgorithm(G3_Runtime*, int, G3_Char**);

std::unordered_map<std::string, EquiSolnAlgo* (*)(G3_Runtime*, int, G3_Char**)> 
EquiSolnAlgoLibrary = {
  {"Linear",            G3Parse_newLinearAlgorithm},
  {"SecantNewton",      G3Parse_newSecantNewtonAlgorithm},
  {"Newton",            G3Parse_newEquiSolnAlgo},
  {"NewtonHallM",       G3Parse_newEquiSolnAlgo},
  {"NewtonHall",        G3Parse_newEquiSolnAlgo},
  {"ModifiedNewton",    G3Parse_newEquiSolnAlgo},
  {"KrylovNewton",      G3Parse_newEquiSolnAlgo},
  {"RaphsonNewton",     G3Parse_newEquiSolnAlgo},
  {"MillerNewton",      G3Parse_newEquiSolnAlgo},
  {"PeriodicNewton",    G3Parse_newEquiSolnAlgo},
  {"Broyden",           G3Parse_newEquiSolnAlgo},
  {"BFGS",              G3Parse_newEquiSolnAlgo},
  {"NewtonLineSearch",  G3Parse_newEquiSolnAlgo},
  {"ExpressNewton",     G3Parse_newEquiSolnAlgo},
};

