#pragma once
#include <stdio.h>
#include <vector>
#include <array>

// ============================================
// VectorND and MatrixND are not required here
// because the linear algebra operations are not used
// --------------------------------------------
#include <VectorND.h>
// #include <MatrixND.h>
using OpenSees::VectorND;
// using OpenSees::MatrixND;
// ============================================
template <int ndim>
class PeriParticle
{
public:
//  using Array = std::array<double, ndim>;
    template <int nd> 
    using Array = VectorND<nd>;

    // // Print a representation of the domain
    // void print(int flag);
    // ============================================
    // DATA
    // ============================================
    int numfam;
    std::vector<int> nodefam;
    std::vector<double> vol;
    std::vector<double> correction, bond_dmg;

    double vol_h;
    Array<ndim> coord;
    Array<ndim> bforce, bdisp, pforce;
    Array<ndim> disp, vel, acc;
    std::array<int, ndim> is_force_bound, is_disp_bound;
    std::array<double, (ndim-1)*3> stress, strain;
    double  energy;
    std::array<double, 6> misc;
};
