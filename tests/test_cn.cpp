// Compares the explicit scheme (bj_flow) with Crank-Nicolson (bj_flow_cn).
// Run from this directory:  make -C .. test   (tables test-*.dat are read from the working dir)
#include "../bj_flow.h"
#include "../bj_flow_cn.h"
#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

static double delta(int i, double eta_max, int Nl)
{
	double sigma2 = 2 * .005;
	double eta = -eta_max + 2 * eta_max * i / (Nl - 1);
	return 1. / sqrt(M_PI * sigma2) * exp(-(eta * eta) / sigma2);
}

typedef std::vector<std::vector<double> > Fields;
static Fields snapshot(bj_flow& s) { Fields f(5); f[0] = s.get_B(); f[1] = s.get_Q(); f[2] = s.get_S(); f[3] = s.get_e(); f[4] = s.get_u(); return f; }
static double norm2(const std::vector<double>& a) { double s = 0; for (double x : a) s += x * x; return sqrt(s); }
static double diff2(const std::vector<double>& a, const std::vector<double>& b) { double s = 0; for (size_t i = 0; i < a.size(); i++) s += (a[i] - b[i]) * (a[i] - b[i]); return sqrt(s); }

int main(int argc, char** argv)
{
	const int Nl = 501; const double eta_max = 5, t1 = 1;
	double t_end = argc > 1 ? atof(argv[1]) : 6.0;
	const char* names[5] = { "B", "Q", "S", "E", "U" };
	std::vector<double> dts = { 0.04, 0.02, 0.01, 0.005 };

	for (int src = 0; src < 5; src++)
	{
		Fields IC(5, std::vector<double>(Nl, 0));
		for (int i = 0; i < Nl; i++) IC[src][i] = delta(i, eta_max, Nl) / t1;

		bj_flow ex("test");
		ex.t = t1; ex.initial_conditions(IC);
		auto t0 = std::chrono::steady_clock::now();
		int nsteps = (int)((t_end - t1) / 2e-5);
		for (int k = 0; k < nsteps; k++) ex.step();
		double t_ex = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
		Fields ref = snapshot(ex);
		printf("source %s: explicit dt=2e-5, %d steps, %.1fs, t=%.5f\n", names[src], nsteps, t_ex, ex.t);

		// self-convergence reference: CN with a much smaller step (independent of the explicit scheme's own error)
		bj_flow_cn fine("test", 0.00125);
		fine.t = t1; fine.initial_conditions(IC); fine.evolve_to(ex.t);
		Fields fref = snapshot(fine);

		double prev = 0;
		for (double dt : dts)
		{
			bj_flow_cn cn("test", dt);
			cn.t = t1; cn.initial_conditions(IC);
			t0 = std::chrono::steady_clock::now();
			cn.evolve_to(ex.t);
			double t_cn = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
			Fields f = snapshot(cn);
			double num = 0, den = 0;
			for (int a = 0; a < 5; a++) { double d = diff2(f[a], ref[a]), n = norm2(ref[a]); num += d * d; den += n * n; }
			double err = sqrt(num / den);
			double num2 = 0, den2 = 0;
			for (int a = 0; a < 5; a++) { double d = diff2(f[a], fref[a]), n = norm2(fref[a]); num2 += d * d; den2 += n * n; }
			printf("   CN dt=%-7g rel.diff to explicit = %.3e   to CN(dt=0.00125) = %.3e  (ratio to prev %.2f)  %.3fs\n",
			       dt, err, sqrt(num2 / den2), prev > 0 ? prev / sqrt(num2 / den2) : 0.0, t_cn);
			err = sqrt(num2 / den2);
			prev = err;
		}
	}
	return 0;
}
