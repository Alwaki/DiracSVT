/*
Project:        Shell evolution of the dirac equation
                
Authors:        Alexander Kiessling
                (2022-2023)

Description:    Header file for solver.cpp. Also serves as 
				function documentation.
*/


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

/*
* brief: structs for passing multiple values of different types simultaneously
*/
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

/*
* brief: Solves the coupled dirac equation at a single point 
*/
std::pair <double, double> pointSolve(double x, double F, double G, const containers::parameters& params, 
										double B, double sigmaV0, double sigmaR, double sigmaa, double dV0, 
										double dR, double da);

/*
* brief: Calls pointSolve repeatedly to integrate along many points
*/
std::pair <double, double> IntegrateRK4(double iniF, double iniG, const containers::parameters& params, 
											double B, double a0, double sigmaV0, double sigmaR, double sigmaa, 
											double dV0, double dR, double da, double xend, double xstart);

/*
* brief: Same as IntegrateRK4, but returns solutions for all solved points
*/
containers::solutionsFG IntegrateRK4_withInfo(double iniF, double iniG, const containers::parameters& params, 
											double B, double a0, double sigmaV0, double sigmaR, double sigmaa, 
											double dV0, double dR, double da, double xend, double xstart);

/*
* brief: Boundary conditions
*/
std::tuple<double, double, double, double> BC(const containers::parameters& params, double B, double a0, double sigmaV0,
											double sigmaR, double sigmaa, double dV0, double dR, double da);

/*
* brief: Positive boundary conditions
*/
std::tuple<double, double, double, double> BC_pos(const containers::parameters& params, double B, double a0, double sigmaV0,
											double sigmaR, double sigmaa, double dV0, double dR, double da);

/*
* brief: Main solver loop to iterate solution for convergence.
*
* note: Returns B=100, a0 = 0 upon divergence.
*/
containers::solutionsBa0 solveDirac(const containers::parameters& params, double a0_in, double B0);

#endif