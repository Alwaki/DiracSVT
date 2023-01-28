#ifndef SOLVER
#define SOLVER

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <tuple>
#include <iostream>
#include <fstream>
#include <limits>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error

#define NEUTRON_MASS_MEV 939.56542052
#define PROTON_MASS_MEV 938.27208816

namespace containers
{
	typedef struct parameters
	{
    public:
        double V0, kappa, lambda, r0, a, Rls, als, N, Z, l, k, xmin, xmax, xmatch,
                m, tensorV, kappa_so;
        int scenario, isospin;

	} parameters;

	typedef struct solutionsFG
	{
		public: 
			std::pair <double, double> FG;
			std::vector<double> rvals;
			std::vector<double> fvals;
			std::vector<double> gvals;
	} solutionsFG;

	typedef struct solutionsBa0
	{
		public: 
			double B;
			double a0;
			std::vector<double> rvals;
			std::vector<double> fvals;
			std::vector<double> gvals;
	} solutionsBa0;
}

std::pair <double, double> pointSolve(double x, double F, double G, const containers::parameters& params, 
										double B, double sigmaV0, double sigmaR, double sigmaa, double dV0, 
										double dR, double da)
{
	double Col = 0;

	if (params.isospin == 1)
	{
		if (x > sigmaR)
		{
			Col = 0.0072923 * params.Z / x;
		}
		else
		{
			Col = 0.0072923 * params.Z * (3 * pow(sigmaR, 2) - pow(x, 2)) / (2 * pow(sigmaR, 3));
		}
	}
	
	double sigma = sigmaV0 / (1 + exp((x - sigmaR) / sigmaa)) + Col;
	double delta = dV0 / (1 + exp((x - dR) / da)) + Col;
	double U = params.tensorV / (1 + exp((x - sigmaR) / sigmaa));
	double dfgW_1 = (-B + sigma) * G + (params.k / x - U) * F;
	double dfgW_2 = (2 * params.m + B - delta) * F + (U - params.k / x) * G;
	return { dfgW_1, dfgW_2 };
}

std::pair <double, double> IntegrateRK4(double iniF, double iniG, const containers::parameters& params, 
											double B, double a0, double sigmaV0, double sigmaR, double sigmaa, 
											double dV0, double dR, double da, double xend, double xstart)
{
	double h = 0.001;
	double step = (xend - xstart) * h;
	double F = iniF;
	double G = iniG;
	double step_2 = step / 2;
	double step_6 = step / 6;
	int end_condition = int(1.0 / h);

	double x = 0;

	std::pair <double, double> v1, v2, v3, v4;

	for (int i = 0; i < end_condition; i++)
	{
		x = i * step + xstart;
		v1 = pointSolve(x, F, G, params, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		v2 = pointSolve(x + step_2, F + v1.first * step_2, G + v1.second * step_2, 
			 params, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		v3 = pointSolve(x + step_2, F + v2.first * step_2, G + v2.second * step_2, 
			 params, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		v4 = pointSolve(x + step, F + v3.first * step, G + v3.second * step, 
			params, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		F += (v1.first + 2 * v2.first + 2 * v3.first + v4.first) * step_6;
		G += (v1.second + 2 * v2.second + 2 * v3.second + v4.second) * step_6;
	}

	return { F, G };
}

containers::solutionsFG IntegrateRK4_withInfo(double iniF, double iniG, const containers::parameters& params, 
											double B, double a0, double sigmaV0, double sigmaR, double sigmaa, 
											double dV0, double dR, double da, double xend, double xstart)
{
	containers::solutionsFG solution;
	
	double h = 0.001;
	double step = (xend - xstart) * h;
	double F = iniF;
	double G = iniG;
	double step_2 = step / 2;
	double step_6 = step / 6;
	int end_condition = int(1.0 / h);

	std::vector<double> rvals;
	std::vector<double> fvals;
	std::vector<double> gvals;

	double x = 0;

	std::pair <double, double> v1, v2, v3, v4;

	for (int i = 0; i < end_condition; i++)
	{
		x = i * step + xstart;
		v1 = pointSolve(x, F, G, params, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		v2 = pointSolve(x + step_2, F + v1.first * step_2, G + v1.second * step_2, 
			 params, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		v3 = pointSolve(x + step_2, F + v2.first * step_2, G + v2.second * step_2, 
			 params, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		v4 = pointSolve(x + step, F + v3.first * step, G + v3.second * step, 
			params, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		F += (v1.first + 2 * v2.first + 2 * v3.first + v4.first) * step_6;
		G += (v1.second + 2 * v2.second + 2 * v3.second + v4.second) * step_6;
		solution.rvals.push_back(x);
		solution.fvals.push_back(F);
		solution.gvals.push_back(G);
	}
	solution.FG = { F, G };

	return solution;
}

std::tuple<double, double, double, double> BC(const containers::parameters& params, double B, double a0, double sigmaV0,
											double sigmaR, double sigmaa, double dV0, double dR, double da)
{
	double Foutbc, Goutbc, miu, Finbc, Ginbc;
	if (params.k < 0)
	{
		Foutbc = -a0 * pow(params.xmin, (params.l + 2)) * (-B + sigmaV0 / (1 + exp(-sigmaR / sigmaa))) / (params.l + 2 - params.k);
		Goutbc = a0 * pow(params.xmin, (params.l + 1));
	}
	else
	{
		Foutbc = a0 * pow(params.xmin, params.l);
		Goutbc = a0 * pow(params.xmin, (params.l + 1)) * (2 * params.m + B - dV0 / (1 + exp(-dR / da))) / (params.l + params.k + 1);

	}
	try
	{
		miu = sqrt(-2 * params.m * B - pow(B, 2));
		Finbc = -sqrt(-B / (2 * params.m + B)) * exp(-miu * params.xmax);
	}
	catch (...)
	{
		miu = 1e-10;
		Finbc = -sqrt(-B / (1e-10)) * exp(-miu * params.xmax);
	}
	Ginbc = exp(-miu * params.xmax);

	return { Foutbc, Goutbc, Finbc, Ginbc };
}

std::tuple<double, double, double, double> BC_pos(const containers::parameters& params, double B, double a0, double sigmaV0,
											double sigmaR, double sigmaa, double dV0, double dR, double da)
{
	double Foutbc, Goutbc, Finbc, Ginbc;
	if (params.k < 0)
	{
		Foutbc = -a0 * pow(params.xmin, (params.l + 2)) * (-B + sigmaV0 / (1 + exp(-sigmaR / sigmaa))) / (params.l + 2 - params.k);
		Goutbc = a0 * pow(params.xmin, (params.l + 1));
		Finbc = (std::sph_bessel(params.l, params.xmax) + std::sph_neumann(params.l, params.xmax));
		Ginbc = sqrt(pow(B, 2) + 2 * B * params.m) / (B + 2 * params.m) * (std::sph_bessel(params.l + 1, params.xmax) + std::sph_neumann(params.l + 1, params.xmax));
	}
	else
	{
		Foutbc = a0 * pow(params.xmin, params.l);
		Goutbc = a0 * pow(params.xmin, (params.l + 1)) * (2 * params.m + B - dV0 / (1 + exp(-dR / da))) / (params.l + params.k + 1);
		Finbc = (std::sph_bessel(params.l, params.xmax) + std::sph_neumann(params.l, params.xmax));
		Ginbc = sqrt(pow(B, 2) + 2 * B * params.m) / (B + 2 * params.m) * (std::sph_bessel(params.l - 1, params.xmax) + std::sph_neumann(params.l - 1, params.xmax));
	}
	double norm = sqrt(pow(std::sph_bessel(params.l, params.xmax), 2) + pow(std::sph_neumann(params.l, params.xmax), 2));

	return { Foutbc, Goutbc, Finbc / norm, Ginbc / norm };
}

containers::solutionsBa0 solveDirac(const containers::parameters& params, double a0_in, double B0)
{
	// Setup potential
	double A = params.N + params.Z;
	double sigmaV0 = 0;
	double dV0 = 0;
	if (params.isospin == 1)
	{
		sigmaV0 = params.V0*(1 + params.kappa*(params.N-params.Z)/A);
		if (params.scenario == 1)
		{
			dV0 = -params.lambda*sigmaV0;
		}
		else if (params.scenario==2)
		{
			dV0 = -params.lambda*params.V0*(1 - params.kappa*(params.N-params.Z)/A);
		}
		else if (params.scenario==3)
		{
			dV0 = -params.lambda*params.V0*(1 - params.kappa_so*(params.N-params.Z)/A);
		}
	}
	else if (params.isospin == -1)
	{
		sigmaV0 = params.V0*(1 - params.kappa*(params.N-params.Z)/A);
		if (params.scenario==1)
		{
			dV0 = -params.lambda * sigmaV0;
		}
		else if (params.scenario==2)
		{
			dV0 = -params.lambda*params.V0*(1 + params.kappa*(params.N-params.Z)/A);
		}
		else if (params.scenario==3)
		{
			dV0 = -params.lambda*params.V0*(1 + params.kappa_so*(params.N-params.Z)/A);
		}
		
		
	}
    double sigmaR = params.r0*pow(A,1.0/3.0);
    double dR = params.Rls*pow(A,1.0/3.0);
    double sigmaa = params.a;
    double da = params.a;
    double error = 100;
	double B = B0;
	double a0 = a0_in;

	// Setup solution container structs
	containers::solutionsFG inSol;
	containers::solutionsFG outSol;
	containers::solutionsBa0 convergedSol;

	// Iterate solvers
	double h = 0.0001;
	int iterations = 0;
	while (error > 0.0001)
	{
		double Foutbc, Goutbc, Finbc, Ginbc;
		if(B <0)
		{
			std::tie(Foutbc, Goutbc, Finbc, Ginbc) = BC(params, B, a0, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		}
		else
		{
			std::tie(Foutbc, Goutbc, Finbc, Ginbc) = BC_pos(params, B, a0, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		}
		inSol = IntegrateRK4_withInfo(Finbc, Ginbc, params, B, a0, sigmaV0, 
										sigmaR, sigmaa, dV0, dR, da, params.xmatch, params.xmax);
        outSol = IntegrateRK4_withInfo(Foutbc, Goutbc, params, B, a0, sigmaV0, 
										sigmaR, sigmaa, dV0, dR, da, params.xmatch, params.xmin);
		double B1 = B + B*h;
		
		std::pair <double, double> inFG = inSol.FG;
		std::pair <double, double> outFG = outSol.FG;

		if(B1 < 0)
		{
			std::tie(Foutbc, Goutbc, Finbc, Ginbc) = BC(params, B1, a0, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		}
		else
		{
			std::tie(Foutbc, Goutbc, Finbc, Ginbc) = BC_pos(params, B1, a0, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		}
		std::pair <double, double> dBinFG = IntegrateRK4(Finbc, Ginbc, params, B1, a0, sigmaV0, 
										sigmaR, sigmaa, dV0, dR, da, params.xmatch, params.xmax);
        std::pair <double, double> dBoutFG = IntegrateRK4(Foutbc, Goutbc, params, B1, a0, sigmaV0, 
										sigmaR, sigmaa, dV0, dR, da, params.xmatch, params.xmin);
										
		double dGFdB_1    = ((dBoutFG.first - dBinFG.first) - (outFG.first - inFG.first))/(B*h);
		double dGFdB_2    = ((dBoutFG.second - dBinFG.second) - (outFG.second - inFG.second))/(B*h);
		double da0outGF_1 = outFG.first * (1.0 + h);
		double da0outGF_2 = outFG.second * (1.0 + h);
		double dGFda0_1   = ((da0outGF_1 - inFG.first) - (outFG.first - inFG.first)) / (a0*h);
		double dGFda0_2   = ((da0outGF_2 - inFG.second) - (outFG.second - inFG.second)) / (a0*h);
		double dOutIn_1   = (outFG.first - inFG.first);
		double dOutIn_2   = (outFG.second - inFG.second);

		Eigen::MatrixXd M(2,2);
		M(0,0) = dGFdB_1;
		M(0,1) = dGFda0_1;
		M(1,0) = dGFdB_2;
		M(1,1) = dGFda0_2;
		Eigen::VectorXd Old(2), Cold(2), diff(2), New(2);
		Old(0) = B;
		Old(1) = a0;
		diff(0) = dOutIn_1, 
		diff(1) = dOutIn_2;
		Cold = (M*Old)-diff;
		New = M.completeOrthogonalDecomposition().solve(Cold);
		B = New(0);
		a0 = New(1);
		error = (New-Old).norm();

		iterations++;
	}
convergedSol.B = B;
convergedSol.a0 = a0;

for(size_t i = 0; i < outSol.rvals.size(); i++)
{
	convergedSol.rvals.push_back(outSol.rvals[i]);
	convergedSol.fvals.push_back(outSol.fvals[i]);
	convergedSol.gvals.push_back(outSol.gvals[i]);
}

for(size_t i = 0; i < inSol.rvals.size(); i++)
{
	convergedSol.rvals.push_back(inSol.rvals[inSol.rvals.size()-i]);
	convergedSol.fvals.push_back(inSol.fvals[inSol.rvals.size()-i]);
	convergedSol.gvals.push_back(inSol.gvals[inSol.rvals.size()-i]);
}

return convergedSol;
}

#endif