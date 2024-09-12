#include <vector>
#include <Matrix.h>
#include <Vector.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <PeriParticle.h>

using OpenSees::VectorND;
using OpenSees::MatrixND;

class PeriDomain {
public:
  //
  // MEMBER FUNCTIONS
  //

  // the special "constructor" member function; no return type declared
  PeriDomain(int totnode, int maxfam);



  // A normal member function
  int hello(double);

  // Print a representation of the domain
  void print(int flag);


  //
  // MEMBER DATA
  //

  // container of particles
  std::vector<PeriParticle<3>> particles;

  int     ndim, totnode, maxfam;
  double  space, delta, vol;
  int     n_bforce, n_bdisp;

  std::vector<int>   numfam;

  std::vector<std::vector<int>> 
      is_force_bound, is_disp_bound, nodefam; // integer matrices

  Matrix    correction, bond_dmg;
  Matrix    bforce_cond, bdisp_cond;
  Matrix    pforce, stress, strain;
  Matrix    misc;

  char plane_type;
};

