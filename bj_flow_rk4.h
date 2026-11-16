#pragma once
#include "bj_flow_cn.h"

// Same linearised Bjorken equations and the same spatial operator A(t) as bj_flow / bj_flow_cn
// (central differences, periodic eta grid), advanced with the classical explicit 4th-order
// Runge-Kutta scheme:   y' = A(t) y,   stages at t, t+dt/2, t+dt/2, t+dt.
//
// RK4 is explicit, so dt is limited by diffusion stability: dt < ~2.78 * h^2 t^2 / (4 max(D, 4 etas/(3 M)))
// (about 9e-4 at t = 1 for the test tables and growing like t^2). Much larger than the 2e-5 of the
// forward-Euler code, and the error is O(dt^4).
class bj_flow_rk4 : public bj_flow_cn
{
public:
	explicit bj_flow_rk4(const std::string& name = "1n0", double dt_rk4_ = 2e-4, bool fix_dQP_in_pressure = false)
		: bj_flow_cn(name, dt_rk4_, fix_dQP_in_pressure) {}

	void step(double dt_step);
	void step() { step(dt_cn); }          // dt_cn doubles as the RK4 step size
	void evolve_to(double t_end);
};
