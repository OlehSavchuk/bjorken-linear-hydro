#include "bj_flow_cn.h"
#include <array>
#include <cmath>
#include <cstring>
#include <algorithm>

namespace
{
typedef double M5[5][5];
typedef double V5[5];

void zero(M5 a) { std::memset(a, 0, sizeof(M5)); }
void copy(M5 dst, const M5 src) { std::memcpy(dst, src, sizeof(M5)); }
void mul(M5 c, const M5 a, const M5 b)           // c = a*b  (c must not alias a or b)
{
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
		{
			double s = 0;
			for (int k = 0; k < 5; k++) s += a[i][k] * b[k][j];
			c[i][j] = s;
		}
}
void matvec(V5 y, const M5 a, const V5 x)        // y = a*x (y must not alias x)
{
	for (int i = 0; i < 5; i++)
	{
		double s = 0;
		for (int k = 0; k < 5; k++) s += a[i][k] * x[k];
		y[i] = s;
	}
}
// b = a^{-1} b for an n-column right-hand side (columns of M5), Gauss-Jordan with partial pivoting
void solve_inplace(const M5 a_in, M5 b)
{
	M5 a;
	copy(a, a_in);
	for (int c = 0; c < 5; c++)
	{
		int p = c;
		for (int r = c + 1; r < 5; r++)
			if (std::fabs(a[r][c]) > std::fabs(a[p][c])) p = r;
		if (p != c)
			for (int j = 0; j < 5; j++) { std::swap(a[c][j], a[p][j]); std::swap(b[c][j], b[p][j]); }
		double inv = 1.0 / a[c][c];
		for (int j = 0; j < 5; j++) { a[c][j] *= inv; b[c][j] *= inv; }
		for (int r = 0; r < 5; r++)
		{
			if (r == c) continue;
			double f = a[r][c];
			if (f == 0) continue;
			for (int j = 0; j < 5; j++) { a[r][j] -= f * a[c][j]; b[r][j] -= f * b[c][j]; }
		}
	}
}
void solve_vec(const M5 a, V5 b)                 // b = a^{-1} b
{
	M5 rhs;
	zero(rhs);
	for (int i = 0; i < 5; i++) rhs[i][0] = b[i];
	solve_inplace(a, rhs);
	for (int i = 0; i < 5; i++) b[i] = rhs[i][0];
}
}

// Field order: 0=B 1=Q 2=S 3=e 4=u.   A y|_i = Lo y_{i-1} + Dg y_i + Up y_{i+1}
void bj_flow_cn::blocks(double t, Mat Lo, Mat Dg, Mat Up)
{
	zero(Lo); zero(Dg); zero(Up);
	const double h = deta, h2 = deta * deta, t2 = t * t;
	const double dBP_ = dBP(t), dQP_ = dQP(t), dSP_ = dSP(t), c2 = cs2(t);
	const double M = Mst(t), dM = dMst(t), eta = etas(t), Dd = D(t);
	const double aB = (rhoB(t) + Dd * dt_rhoB(t)) / t;
	const double aQ = (rhoQ(t) + Dd * dt_rhoQ(t)) / t;
	const double bE = (M - 4. * eta / (3. * t)) / t;
	const double d = Dd / t2;                       // diffusion coefficient of B, Q, S
	const double c = 4. / 3. * eta / t2 / M;       // viscous coefficient in the u equation

	for (int f = 0; f < 3; f++)                     // B, Q, S: -f/t + d f''
	{
		Dg[f][f] = -1. / t - 2. * d / h2;
		Lo[f][f] = Up[f][f] = d / h2;
	}
	Lo[0][4] = -aB / (2 * h);  Up[0][4] = aB / (2 * h);     // + aB * u'
	Lo[1][4] = -aQ / (2 * h);  Up[1][4] = aQ / (2 * h);

	// e:  -(1 + cs2) e / t - (dBP B + dP_Q Q + dSP S)/t - bE u'
	Dg[3][3] = -(1. + c2) / t;
	Dg[3][0] = -dBP_ / t;
	Dg[3][1] = -(fix_dQP ? dQP_ : dBP_) / t;
	Dg[3][2] = -dSP_ / t;
	Lo[3][4] = bE / (2 * h);   Up[3][4] = -bE / (2 * h);

	// u:  -(2/t + M'/M) u - (cs2 e' + dBP B' + dQP Q' + dSP S')/(M t) + c u''
	Dg[4][4] = -(2. / t + dM / M) - 2. * c / h2;
	Lo[4][4] = Up[4][4] = c / h2;
	const double g[4] = { dBP_, dQP_, dSP_, c2 };            // coefficients of B, Q, S, e gradients
	const int col[4] = { 0, 1, 2, 3 };
	for (int k = 0; k < 4; k++)
	{
		double w = g[k] / (M * t) / (2 * h);
		Lo[4][col[k]] = +w;
		Up[4][col[k]] = -w;
	}
}

void bj_flow_cn::step(double dts)
{
	const int N = Nl;
	const double tm = t + 0.5 * dts;
	M5 Lo, Dg, Up;
	blocks(tm, Lo, Dg, Up);
	const double hdt = 0.5 * dts;

	// current state as 5-vectors
	std::vector<std::array<double, 5> > y(N), r(N);
	for (int i = 0; i < N; i++) { y[i][0] = B[i]; y[i][1] = Q[i]; y[i][2] = S[i]; y[i][3] = e[i]; y[i][4] = u[i]; }

	// r = (I + dt/2 A) y
	for (int i = 0; i < N; i++)
	{
		const std::array<double, 5>& ym = y[(i + N - 1) % N];
		const std::array<double, 5>& yp = y[(i + 1) % N];
		for (int a = 0; a < 5; a++)
		{
			double s = y[i][a];
			for (int b = 0; b < 5; b++)
				s += hdt * (Lo[a][b] * ym[b] + Dg[a][b] * y[i][b] + Up[a][b] * yp[b]);
			r[i][a] = s;
		}
	}

	// implicit matrix blocks:  Lm y_{i-1} + Dm y_i + Um y_{i+1} = r_i
	M5 Lm, Dm, Um;
	for (int a = 0; a < 5; a++)
		for (int b = 0; b < 5; b++)
		{
			Lm[a][b] = -hdt * Lo[a][b];
			Um[a][b] = -hdt * Up[a][b];
			Dm[a][b] = (a == b ? 1. : 0.) - hdt * Dg[a][b];
		}

	// periodic block-Thomas: y_i = a_i + P_i w + Q_i y_{i+1}, with w = y_{N-1} unknown.
	std::vector<std::array<double, 5> > av(N), sv(N);
	std::vector<std::array<double, 25> > P(N), Qm(N), T(N);
	auto asM = [](std::array<double, 25>& arr) { return reinterpret_cast<M5*>(arr.data()); };

	M5 G, tmp, tmp2;
	V5 vtmp, vtmp2;
	for (int i = 0; i < N - 1; i++)
	{
		// G = Dm + Lm Q_{i-1}
		copy(G, Dm);
		if (i > 0)
		{
			mul(tmp, Lm, *asM(Qm[i - 1]));
			for (int a = 0; a < 5; a++) for (int b = 0; b < 5; b++) G[a][b] += tmp[a][b];
		}
		// a_i = G^{-1}(r_i - Lm a_{i-1})
		for (int a = 0; a < 5; a++) vtmp[a] = r[i][a];
		if (i > 0)
		{
			V5 prev; for (int a = 0; a < 5; a++) prev[a] = av[i - 1][a];
			matvec(vtmp2, Lm, prev);
			for (int a = 0; a < 5; a++) vtmp[a] -= vtmp2[a];
		}
		solve_vec(G, vtmp);
		for (int a = 0; a < 5; a++) av[i][a] = vtmp[a];
		// Q_i = -G^{-1} Um
		for (int a = 0; a < 5; a++) for (int b = 0; b < 5; b++) tmp[a][b] = -Um[a][b];
		solve_inplace(G, tmp);
		std::memcpy(Qm[i].data(), tmp, sizeof(M5));
		// P_i = -G^{-1} Lm P_{i-1}   (P_0 = -G^{-1} Lm)
		if (i == 0) { for (int a = 0; a < 5; a++) for (int b = 0; b < 5; b++) tmp[a][b] = -Lm[a][b]; }
		else { mul(tmp, Lm, *asM(P[i - 1])); for (int a = 0; a < 5; a++) for (int b = 0; b < 5; b++) tmp[a][b] = -tmp[a][b]; }
		solve_inplace(G, tmp);
		std::memcpy(P[i].data(), tmp, sizeof(M5));
	}
	// back-substitution: y_i = s_i + T_i w, i = N-2 ... 0
	{
		int i = N - 2;
		sv[i] = av[i];
		for (int a = 0; a < 5; a++) for (int b = 0; b < 5; b++) tmp[a][b] = (*asM(P[i]))[a][b] + (*asM(Qm[i]))[a][b];
		std::memcpy(T[i].data(), tmp, sizeof(M5));
		for (i = N - 3; i >= 0; i--)
		{
			V5 snext; for (int a = 0; a < 5; a++) snext[a] = sv[i + 1][a];
			matvec(vtmp, *asM(Qm[i]), snext);
			for (int a = 0; a < 5; a++) sv[i][a] = av[i][a] + vtmp[a];
			mul(tmp, *asM(Qm[i]), *asM(T[i + 1]));
			for (int a = 0; a < 5; a++) for (int b = 0; b < 5; b++) tmp[a][b] += (*asM(P[i]))[a][b];
			std::memcpy(T[i].data(), tmp, sizeof(M5));
		}
	}
	// closing row N-1:  Lm y_{N-2} + Dm w + Um y_0 = r_{N-1}
	{
		M5 Cm;
		mul(tmp, Lm, *asM(T[N - 2]));
		mul(tmp2, Um, *asM(T[0]));
		for (int a = 0; a < 5; a++) for (int b = 0; b < 5; b++) Cm[a][b] = Dm[a][b] + tmp[a][b] + tmp2[a][b];
		V5 rhs, s1, s2;
		for (int a = 0; a < 5; a++) { rhs[a] = r[N - 1][a]; s1[a] = sv[N - 2][a]; s2[a] = sv[0][a]; }
		matvec(vtmp, Lm, s1);
		matvec(vtmp2, Um, s2);
		for (int a = 0; a < 5; a++) rhs[a] -= vtmp[a] + vtmp2[a];
		solve_vec(Cm, rhs);                        // rhs := w
		for (int i = 0; i < N - 1; i++)
		{
			V5 w; for (int a = 0; a < 5; a++) w[a] = rhs[a];
			matvec(vtmp, *asM(T[i]), w);
			for (int a = 0; a < 5; a++) y[i][a] = sv[i][a] + vtmp[a];
		}
		for (int a = 0; a < 5; a++) y[N - 1][a] = rhs[a];
	}

	for (int i = 0; i < N; i++) { B[i] = y[i][0]; Q[i] = y[i][1]; S[i] = y[i][2]; e[i] = y[i][3]; u[i] = y[i][4]; }
	t += dts;
}

void bj_flow_cn::evolve_to(double t_end)
{
	while (t_end - t > 1e-12)
		step(std::min(dt_cn, t_end - t));
}
