#pragma once
#include <stdio.h>
#include <VectorND.h>
#include <MatrixND.h>

using OpenSees::VectorND;
using OpenSees::MatrixND;

template<int ndm>
class PeriParticle {
public:

  // Print a representation of the domain
  void print(int flag);

  // DATA
  VectorND<ndm> coord, 
                pforce, bforce, bdisp, 
                disp, vel, acc;

  double energy, vol_h, vol;

};

//
// Implementation
//
template <int ndm>
void
PeriParticle<ndm>::print(int flag)
{

  printf("{\"Coordinates\": [");
  for (int i=0; i<ndm; i++)
    printf("%f ", coord[i]);
  printf("]}\n");
}
