#pragma once
#include <stdio.h>
#include <vector>
#include <array>

// ============================================
// VectorND and MatrixND are not required here
// because the linear algebra operations are not used
// --------------------------------------------
// #include <VectorND.h>
// #include <MatrixND.h>
// using OpenSees::VectorND;
// using OpenSees::MatrixND;
// ============================================
template <int ndim>
class PeriParticle
{
public:
	// // Print a representation of the domain
	// void print(int flag);
	// ============================================
	// DATA
	// ============================================
	int numfam;
	std::vector<int> nodefam;
	std::array<double, ndim> coord;
	std::vector<double> vol;
	double vol_h;
	std::vector<double> correction, bond_dmg;
	std::array<int, ndim> is_force_bound, is_disp_bound;
	std::array<double, ndim> bforce, bdisp, pforce;
	std::array<double, (ndim - 1) * 3> stress, strain;
	std::array<double, ndim> disp, vel, acc;
	double energy;
	std::array<double, 6> misc;
};
