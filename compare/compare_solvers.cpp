// Green's functions G_{A->B}(eta) at t_end for a delta source in field A at start time t1, computed with
//   euler : forward Euler (bj_flow, dt = 2e-5, the original code)
//   cn    : Crank-Nicolson (dt = 0.01)
//   rk4   : classical Runge-Kutta 4 (dt = 2e-4)
// Output (in out_dir):  green_<solver>_t1_<t1>_src<A>.dat   5 rows (responses B,Q,S,E,Uz) x Nl columns
//                       params.txt                          run parameters (read by plot_green_comparison.py)
//   usage:  compare_solvers <table_prefix> <out_dir>      (run from the folder holding the <prefix>-*.dat tables)
#include "../bj_flow.h"
#include "../bj_flow_cn.h"
#include "../bj_flow_rk4.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

static const int Nl = 501;
static const double eta_max = 5., t_end = 11., PI_DRIVER = 3.14;   // PI = 3.14 as in Hydro-fluks.cpp

static double delta(int i)
{
	double sigma2 = 2 * .005, eta = -eta_max + 2 * eta_max * i / (Nl - 1);
	return 1. / sqrt(PI_DRIVER * sigma2) * exp(-(eta * eta) / sigma2);
}

typedef std::vector<std::vector<double> > Fields;

static void write_rows(const std::string& path, bj_flow& s)
{
	std::ofstream f(path);
	std::vector<double> rows[5] = { s.get_B(), s.get_Q(), s.get_S(), s.get_e(), s.get_u() };
	f.precision(12);
	for (int r = 0; r < 5; r++) { for (int j = 0; j < Nl; j++) f << rows[r][j] << (j + 1 < Nl ? "\t" : ""); f << "\n"; }
}

int main(int argc, char** argv)
{
	std::string prefix = argc > 1 ? argv[1] : "1n0", out = argc > 2 ? argv[2] : ".";
	const char* src_name[5] = { "B", "Q", "S", "E", "U" };
	const double t1s[2] = { 1.0, 6.0 };
	const double dt_euler = 2e-5, dt_cn = 0.01, dt_rk4 = 5e-4;
	std::map<std::string, double> seconds;

	for (double t1 : t1s)
		for (int src = 0; src < 5; src++)
		{
			Fields IC(5, std::vector<double>(Nl, 0));
			for (int i = 0; i < Nl; i++) IC[src][i] = delta(i) / t1;
			char tag[64]; snprintf(tag, sizeof tag, "t1_%g_src%s", t1, src_name[src]);

			auto timed = [&](const char* name, bj_flow& s, auto&& run)
			{
				auto t0 = std::chrono::steady_clock::now();
				run();
				seconds[name] += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
				write_rows(out + "/green_" + name + "_" + tag + ".dat", s);
			};

			bj_flow ex(prefix); ex.t = t1; ex.initial_conditions(IC);
			timed("euler", ex, [&] { int n = (int)llround((t_end - t1) / dt_euler); for (int k = 0; k < n; k++) ex.step(); });
			bj_flow_cn cn(prefix, dt_cn); cn.t = t1; cn.initial_conditions(IC);
			timed("cn", cn, [&] { cn.evolve_to(t_end); });
			bj_flow_rk4 rk(prefix, dt_rk4); rk.t = t1; rk.initial_conditions(IC);
			timed("rk4", rk, [&] { rk.evolve_to(t_end); });
			fprintf(stderr, "done %s\n", tag);
		}

	std::ofstream p(out + "/params.txt");
	p << "table_prefix=" << prefix << "\nNl=" << Nl << "\neta_max=" << eta_max << "\nt_end=" << t_end
	  << "\nt1_values=1,6\nsource_sigma2=0.01\nsource_PI=" << PI_DRIVER
	  << "\ndt_euler=" << dt_euler << "\ndt_cn=" << dt_cn << "\ndt_rk4=" << dt_rk4
	  << "\nseconds_euler=" << seconds["euler"] << "\nseconds_cn=" << seconds["cn"] << "\nseconds_rk4=" << seconds["rk4"] << "\n";
	return 0;
}
