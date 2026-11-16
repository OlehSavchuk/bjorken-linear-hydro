// Green's functions of the linearised Bjorken hydro with the Crank-Nicolson solver.
// Same output as Hydro-fluks.cpp (files Green-<src>-<resp>.dat, one line per start time, Nl values per line),
// but ~10x faster and with a user-chosen time step.
//
//   ./Hydro-fluks-cn [dt=0.01] [Ntstart=500] [table_prefix=1n0] [out_dir=.]
#include "bj_flow_cn.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

static double delta(int i, double eta_max, int Nl)
{
	const double PI = 3.14;           // as in Hydro-fluks.cpp, so that the two programs agree; M_PI is more exact
	double sigma2 = 2 * .005;
	double eta = -eta_max + 2 * eta_max * i / (Nl - 1);
	return 1. / sqrt(PI * sigma2) * exp(-(eta * eta) / sigma2);
}

int main(int argc, char** argv)
{
	double dt = argc > 1 ? atof(argv[1]) : 0.01;
	int Ntstart = argc > 2 ? atoi(argv[2]) : 500;
	std::string prefix = argc > 3 ? argv[3] : "1n0";
	std::string out = argc > 4 ? argv[4] : ".";
	const int Nl = 501;
	const double eta_max = 5., tstart = 1., tend = 11.;
	const char* src_name[5] = { "B", "Q", "S", "E", "U" };      // IC index 0..4 = B, Q, S, e, u
	const char* resp_name[5] = { "B", "Q", "S", "E", "Uz" };

	bj_flow_cn system(prefix, dt);
	for (int src = 0; src < 5; src++)
	{
		std::vector<std::ofstream> files(5);
		for (int r = 0; r < 5; r++)
			files[r].open(out + "/Green-" + src_name[src] + "-" + resp_name[r] + ".dat");
		for (int ti = 0; ti < Ntstart; ti++)
		{
			double t1 = tstart + (tend - tstart) * ti / (Ntstart - 1);
			system.t = t1;
			std::vector<std::vector<double> > IC(5, std::vector<double>(Nl, 0));
			for (int i = 0; i < Nl; i++) IC[src][i] = delta(i, eta_max, Nl) / t1;
			system.initial_conditions(IC);
			system.evolve_to(tend);
			std::vector<double> res[5] = { system.get_B(), system.get_Q(), system.get_S(), system.get_e(), system.get_u() };
			for (int r = 0; r < 5; r++)
			{
				for (int j = 0; j < Nl; j++) files[r] << res[r][j] << "\t";
				files[r] << std::endl;
			}
		}
		std::cerr << "source " << src_name[src] << " done" << std::endl;
	}
	return 0;
}
