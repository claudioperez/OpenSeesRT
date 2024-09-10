#pragma once
#include <stdio.h>
#include <vector>

// ============================================
// VectorND and MatrixND are not required here
// because the linear algebra operations are not used
// --------------------------------------------
// #include <VectorND.h>
// #include <MatrixND.h>
// using OpenSees::VectorND;
// using OpenSees::MatrixND;
// ============================================

class PeriParticle {
public:

  // // Print a representation of the domain
  // void print(int flag);
  // ============================================
  // DATA
  // ============================================
  int                 numfam;
  std::vector<int>    nodefam;
  std::vector<double> coord, vol;
  double              vol_h;
  std::vector<double> correction, bond_dmg;
  std::vector<int>    is_force_bound, is_disp_bound;
  std::vector<double> bforce, bdisp, pforce, stress, strain,
                      disp, vel, acc;
  double              energy;
  std::vector<double> misc;
};

//
// Implementation
//
// template <int ndm>
// void
// PeriParticle<ndm>::print(int flag)
// {

//   printf("{\"Coordinates\": [");
//   for (int i=0; i<ndm; i++)
//     printf("%f ", coord[i]);
//   printf("]}\n");
// }
