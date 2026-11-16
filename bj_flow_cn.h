#pragma once
#include "bj_flow.h"
#include <array>

// Linearised Bjorken hydro, same equations as bj_flow::step(), but advanced with the
// implicit, second-order Crank-Nicolson scheme
//
//     (I - dt/2 A(t+dt/2)) y^{n+1} = (I + dt/2 A(t+dt/2)) y^n ,   y = (B, Q, S, e, u)
//
// A is the spatial operator with central differences on the periodic eta grid (same stencil
// as the explicit code). Because A has constant coefficients in eta, the implicit system is a
// periodic block-tridiagonal system with 5x5 blocks, solved directly each step in O(Nl).
// The scheme is unconditionally stable, so dt is chosen for accuracy, not stability
// (the explicit code needs dt ~ 2e-5; CN is accurate with dt ~ 1e-2).
class bj_flow_cn : public bj_flow
{
public:
	explicit bj_flow_cn(const std::string& name = "1n0", double dt_cn_ = 0.01, bool fix_dQP_in_pressure = true)
		: bj_flow(name), dt_cn(dt_cn_), fix_dQP(fix_dQP_in_pressure) {}

	double dt_cn;   // time step; the last step is shortened to hit t_end exactly
	// dP = cs2 e + dBP B + dQP Q + dSP S. Older versions of the explicit code used dBP in front of Q (typo);
	// true (default) uses dQP, false reproduces that old behaviour.
	bool fix_dQP;

	void step(double dt_step);   // one CN step of length dt_step
	void step() { step(dt_cn); }
	void evolve_to(double t_end); // advances t -> t_end with steps of at most dt_cn

protected:   // shared with bj_flow_rk4
	typedef double Mat[5][5];
	void blocks(double tm, Mat Lo, Mat Dg, Mat Up);
};
