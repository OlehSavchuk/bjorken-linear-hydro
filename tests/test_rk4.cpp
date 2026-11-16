// RK4 solver: difference to RK4 with a very small step (time-integration error) and to the forward-Euler reference.
// Run from tests/ :  make -C .. test-rk4
#include "../bj_flow.h"
#include "../bj_flow_rk4.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

static double delta(int i, double eta_max, int Nl)
{
	double sigma2 = 2 * .005, eta = -eta_max + 2 * eta_max * i / (Nl - 1);
	return 1. / sqrt(M_PI * sigma2) * exp(-(eta * eta) / sigma2);
}
typedef std::vector<std::vector<double> > Fields;
static Fields snap(bj_flow& s) { Fields f(5); f[0] = s.get_B(); f[1] = s.get_Q(); f[2] = s.get_S(); f[3] = s.get_e(); f[4] = s.get_u(); return f; }
static double rel(const Fields& a, const Fields& b)
{
	double n = 0, d = 0;
	for (int k = 0; k < 5; k++) for (size_t i = 0; i < a[k].size(); i++) { n += (a[k][i] - b[k][i]) * (a[k][i] - b[k][i]); d += b[k][i] * b[k][i]; }
	return sqrt(n / d);
}

int main(int argc, char** argv)
{
	const int Nl = 501; const double eta_max = 5, t1 = 1, t_end = argc > 1 ? atof(argv[1]) : 6.0;
	const char* names[5] = { "B", "Q", "S", "E", "U" };
	for (int src = 0; src < 5; src++)
	{
		Fields IC(5, std::vector<double>(Nl, 0));
		for (int i = 0; i < Nl; i++) IC[src][i] = delta(i, eta_max, Nl) / t1;
		bj_flow ex("test"); ex.t = t1; ex.initial_conditions(IC);
		int n = (int)((t_end - t1) / 2e-5);
		for (int k = 0; k < n; k++) ex.step();
		Fields eul = snap(ex);
		bj_flow_rk4 fine("test", 5e-5); fine.t = t1; fine.initial_conditions(IC); fine.evolve_to(ex.t);
		Fields fref = snap(fine);
		printf("source %s (t = %.3f)\n", names[src], ex.t);
		double prev = 0;
		for (double dt : { 1.1e-3, 8e-4, 4e-4, 2e-4 })
		{
			bj_flow_rk4 s("test", dt); s.t = t1; s.initial_conditions(IC);
			auto t0 = std::chrono::steady_clock::now(); s.evolve_to(ex.t);
			double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
			Fields f = snap(s);
			double e_self = rel(f, fref);
			printf("   RK4 dt=%-7g to RK4(5e-5) %.3e   to Euler(2e-5) %.3e   %.2fs\n", dt, e_self, rel(f, eul), sec);
			prev = e_self;
		}
	}
	return 0;
}
