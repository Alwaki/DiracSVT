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
                m, tensorV, kappa_so, B, a0;
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
										double dR, double da);

std::pair <double, double> IntegrateRK4(double iniF, double iniG, const containers::parameters& params, 
											double B, double a0, double sigmaV0, double sigmaR, double sigmaa, 
											double dV0, double dR, double da, double xend, double xstart);

containers::solutionsFG IntegrateRK4_withInfo(double iniF, double iniG, const containers::parameters& params, 
											double B, double a0, double sigmaV0, double sigmaR, double sigmaa, 
											double dV0, double dR, double da, double xend, double xstart);

std::tuple<double, double, double, double> BC(const containers::parameters& params, double B, double a0, double sigmaV0,
											double sigmaR, double sigmaa, double dV0, double dR, double da);

std::tuple<double, double, double, double> BC_pos(const containers::parameters& params, double B, double a0, double sigmaV0,
											double sigmaR, double sigmaa, double dV0, double dR, double da);

containers::solutionsBa0 solveDirac(const containers::parameters& params, double a0_in, double B0);

#endif