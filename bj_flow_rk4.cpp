#include "bj_flow_rk4.h"
#include <algorithm>
#include <array>
#include <vector>

namespace
{
typedef std::vector<std::array<double, 5> > State;
}

void bj_flow_rk4::step(double dts)
{
	const int N = Nl;
	Mat Lo[3], Dg[3], Up[3];
	const double tt[3] = { t, t + 0.5 * dts, t + dts };       // stage times: t, t+dt/2 (used twice), t+dt
	for (int s = 0; s < 3; s++) blocks(tt[s], Lo[s], Dg[s], Up[s]);

	State y(N), k1(N), k2(N), k3(N), k4(N), ys(N);
	for (int i = 0; i < N; i++) { y[i][0] = B[i]; y[i][1] = Q[i]; y[i][2] = S[i]; y[i][3] = e[i]; y[i][4] = u[i]; }

	auto rhs = [&](int s, const State& in, State& out)
	{
		for (int i = 0; i < N; i++)
		{
			const std::array<double, 5>& ym = in[(i + N - 1) % N];
			const std::array<double, 5>& yp = in[(i + 1) % N];
			for (int a = 0; a < 5; a++)
			{
				double v = 0;
				for (int b = 0; b < 5; b++) v += Lo[s][a][b] * ym[b] + Dg[s][a][b] * in[i][b] + Up[s][a][b] * yp[b];
				out[i][a] = v;
			}
		}
	};
	auto axpy = [&](const State& base, const State& k, double f, State& out)
	{
		for (int i = 0; i < N; i++) for (int a = 0; a < 5; a++) out[i][a] = base[i][a] + f * k[i][a];
	};

	rhs(0, y, k1);
	axpy(y, k1, 0.5 * dts, ys); rhs(1, ys, k2);
	axpy(y, k2, 0.5 * dts, ys); rhs(1, ys, k3);
	axpy(y, k3, dts, ys);       rhs(2, ys, k4);

	for (int i = 0; i < N; i++)
	{
		double v[5];
		for (int a = 0; a < 5; a++) v[a] = y[i][a] + dts / 6. * (k1[i][a] + 2 * k2[i][a] + 2 * k3[i][a] + k4[i][a]);
		B[i] = v[0]; Q[i] = v[1]; S[i] = v[2]; e[i] = v[3]; u[i] = v[4];
	}
	t += dts;
}

void bj_flow_rk4::evolve_to(double t_end)
{
	while (t_end - t > 1e-12)
		step(std::min(dt_cn, t_end - t));
}
