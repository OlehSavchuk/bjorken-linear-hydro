
#include "bj_flow.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>


double delta(int i, double eta_max, int Nl) 
{

	double PI = 3.14;
	double sigma2 = 2*.005;
	double eta = -eta_max + 2 * eta_max * i / (Nl-1);
	return 1. / sqrt(PI*sigma2) * exp(-(eta * eta) / (sigma2));
}

int main()
{
	int Nl = 501;
	double eta_max=5.;
	//for density

	std::ofstream file11;
	std::ofstream file12;
	std::ofstream file13;
	std::ofstream file2;
	std::ofstream file3;
	int Ntstart = 500;
	double tstart = 1;
	double tend = 11;
	bj_flow system;


	file11.open("Green-B-B.dat");
	file12.open("Green-B-Q.dat");
	file13.open("Green-B-S.dat");
	file2.open("Green-B-E.dat");
	file3.open("Green-B-Uz.dat");
	for (int ti = 0; ti < Ntstart; ti++)
	{
		double t1 = tstart +(tend - tstart) * ti / (Ntstart - 1);
		system.t = t1;
		std::vector<std::vector<double>> IC(5, std::vector<double>(Nl, 0));
		for (int i = 0; i < Nl; i++)
		{
			IC[0][i] = delta(i, eta_max, Nl) / t1;
		}
		system.initial_conditions(IC);
		double dt = 10. / 500 * 0.001;
		int Nt = (tend-t1) / dt;
		int Noutput = Nt / 10;
		for (int i = 0; i < (tend - t1) / dt; i++)
		{
			system.step();
		}

			std::vector<double> t_B = system.get_B();
			std::vector<double> t_Q = system.get_Q();
			std::vector<double> t_S = system.get_S();
			std::vector<double> t_e = system.get_e();
			std::vector<double> t_u = system.get_u();
			for (int j = 0; j < Nl; j++)
			{
				file11 << t_B[j] << "\t";
			}
			file11 << std::endl;
			for (int j = 0; j < Nl; j++)
			{
				file12 << t_Q[j] << "\t";
			}
			file12 << std::endl;
			for (int j = 0; j < Nl; j++)
			{
				file13 << t_S[j] << "\t";
			}
			file13 << std::endl;
			for (int j = 0; j < Nl; j++)
			{
				file2 << t_e[j] << "\t";
			}
			file2 << std::endl;
			for (int j = 0; j < Nl; j++)
			{
				file3 << t_u[j] << "\t";
			}
		    file3 << std::endl;

	}


	file11.close();
	file12.close();
	file13.close();
	file2.close();
	file3.close();






	file11.open("Green-Q-B.dat");
	file12.open("Green-Q-Q.dat");
	file13.open("Green-Q-S.dat");
	file2.open("Green-Q-E.dat");
	file3.open("Green-Q-Uz.dat");
	for (int ti = 0; ti < Ntstart; ti++)
	{
		double t1 = tstart + (tend - tstart) * ti / (Ntstart - 1);
		system.t = t1;
		std::vector<std::vector<double>> IC(5, std::vector<double>(Nl, 0));
		for (int i = 0; i < Nl; i++)
		{
			IC[1][i] = delta(i, eta_max, Nl) / t1;
		}
		system.initial_conditions(IC);
		double dt = 10. / 500 * 0.001;
		int Nt = (tend - t1) / dt;
		int Noutput = Nt / 10;
		for (int i = 0; i < (tend - t1) / dt; i++)
		{
			system.step();
		}

		std::vector<double> t_B = system.get_B();
		std::vector<double> t_Q = system.get_Q();
		std::vector<double> t_S = system.get_S();
		std::vector<double> t_e = system.get_e();
		std::vector<double> t_u = system.get_u();
		for (int j = 0; j < Nl; j++)
		{
			file11 << t_B[j] << "\t";
		}
		file11 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file12 << t_Q[j] << "\t";
		}
		file12 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file13 << t_S[j] << "\t";
		}
		file13 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file2 << t_e[j] << "\t";
		}
		file2 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file3 << t_u[j] << "\t";
		}
		file3 << std::endl;

	}


	file11.close();
	file12.close();
	file13.close();
	file2.close();
	file3.close();
		

	file11.open("Green-S-B.dat");
	file12.open("Green-S-Q.dat");
	file13.open("Green-S-S.dat");
	file2.open("Green-S-E.dat");
	file3.open("Green-S-Uz.dat");
	for (int ti = 0; ti < Ntstart; ti++)
	{
		double t1 = tstart + (tend - tstart) * ti / (Ntstart - 1);
		system.t = t1;
		std::vector<std::vector<double>> IC(5, std::vector<double>(Nl, 0));
		for (int i = 0; i < Nl; i++)
		{
			IC[2][i] = delta(i, eta_max, Nl) / t1;
		}
		system.initial_conditions(IC);
		double dt = 10. / 500 * 0.001;
		int Nt = (tend - t1) / dt;
		int Noutput = Nt / 10;
		for (int i = 0; i < (tend - t1) / dt; i++)
		{
			system.step();
		}

		std::vector<double> t_B = system.get_B();
		std::vector<double> t_Q = system.get_Q();
		std::vector<double> t_S = system.get_S();
		std::vector<double> t_e = system.get_e();
		std::vector<double> t_u = system.get_u();
		for (int j = 0; j < Nl; j++)
		{
			file11 << t_B[j] << "\t";
		}
		file11 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file12 << t_Q[j] << "\t";
		}
		file12 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file13 << t_S[j] << "\t";
		}
		file13 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file2 << t_e[j] << "\t";
		}
		file2 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file3 << t_u[j] << "\t";
		}
		file3 << std::endl;

	}


	file11.close();
	file12.close();
	file13.close();
	file2.close();
	file3.close();



	file11.open("Green-E-B.dat");
	file12.open("Green-E-Q.dat");
	file13.open("Green-E-S.dat");
	file2.open("Green-E-E.dat");
	file3.open("Green-E-Uz.dat");
	for (int ti = 0; ti < Ntstart; ti++)
	{
		double t1 = tstart + (tend - tstart) * ti / (Ntstart - 1);
		system.t = t1;
		std::vector<std::vector<double>> IC(5, std::vector<double>(Nl, 0));
		for (int i = 0; i < Nl; i++)
		{
			IC[3][i] = delta(i, eta_max, Nl) / t1;
		}
		system.initial_conditions(IC);
		double dt = 10. / 500 * 0.001;
		int Nt = (tend - t1) / dt;
		int Noutput = Nt / 10;
		for (int i = 0; i < (tend - t1) / dt; i++)
		{
			system.step();
		}

		std::vector<double> t_B = system.get_B();
		std::vector<double> t_Q = system.get_Q();
		std::vector<double> t_S = system.get_S();
		std::vector<double> t_e = system.get_e();
		std::vector<double> t_u = system.get_u();
		for (int j = 0; j < Nl; j++)
		{
			file11 << t_B[j] << "\t";
		}
		file11 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file12 << t_Q[j] << "\t";
		}
		file12 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file13 << t_S[j] << "\t";
		}
		file13 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file2 << t_e[j] << "\t";
		}
		file2 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file3 << t_u[j] << "\t";
		}
		file3 << std::endl;

	}


	file11.close();
	file12.close();
	file13.close();
	file2.close();
	file3.close();



	file11.open("Green-U-B.dat");
	file12.open("Green-U-Q.dat");
	file13.open("Green-U-S.dat");
	file2.open("Green-U-E.dat");
	file3.open("Green-U-Uz.dat");
	for (int ti = 0; ti < Ntstart; ti++)
	{
		double t1 = tstart + (tend - tstart) * ti / (Ntstart - 1);
		system.t = t1;
		std::vector<std::vector<double>> IC(5, std::vector<double>(Nl, 0));
		for (int i = 0; i < Nl; i++)
		{
			IC[4][i] = delta(i, eta_max, Nl) / t1;
		}
		system.initial_conditions(IC);
		double dt = 10. / 500 * 0.001;
		int Nt = (tend - t1) / dt;
		int Noutput = Nt / 10;
		for (int i = 0; i < (tend - t1) / dt; i++)
		{
			system.step();
		}

		std::vector<double> t_B = system.get_B();
		std::vector<double> t_Q = system.get_Q();
		std::vector<double> t_S = system.get_S();
		std::vector<double> t_e = system.get_e();
		std::vector<double> t_u = system.get_u();
		for (int j = 0; j < Nl; j++)
		{
			file11 << t_B[j] << "\t";
		}
		file11 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file12 << t_Q[j] << "\t";
		}
		file12 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file13 << t_S[j] << "\t";
		}
		file13 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file2 << t_e[j] << "\t";
		}
		file2 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file3 << t_u[j] << "\t";
		}
		file3 << std::endl;

	}


	file11.close();
	file12.close();
	file13.close();
	file2.close();
	file3.close();
}