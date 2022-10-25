#include "solver.h"
#include "util.h"

int main()
{
	// Allow user to select scenario and state
	std::cout << "The program requires a scenario and state. These \n";
	std::cout << "are both integers, and the scenario ranges from 1 to 3 \n";
	std::cout << "while the state ranges from 1 to 89. \n";
	std::cout << "\n";
	int scenario = read_user_input("scenario"); 
	int state = read_user_input("state");

	// Load data from files
	auto data = read_csv("data.csv");
	auto parameters = read_csv("parameters.csv");

	// Set parameters dependent on scenario
	double V0 = parameters[scenario][1];
	double kappa = parameters[scenario][2];
	double lambda = parameters[scenario][3];
	double r0 = parameters[scenario][4];
	double a = parameters[scenario][5];
	double Rls = parameters[scenario][6];
	double als = parameters[scenario][7];

	// Set parameters dependent on state
	int isospin = data[state][2];
	double N = data[state][3];
	double Z = data[state][4];
	double l = data[state][5];
	double k = data[state][6];
	double B = data[state][7];
	double xmin = data[state][8];
	double xmax = data[state][9];
	double xmatch = data[state][10];
	double a0 = data[state][11];

	// Set other parameters
	double tensorV = 0;
	double kappa_so = 0;
	double m = 0;
	if (isospin == -1)
		m = NEUTRON_MASS_MEV;
	else
		m = PROTON_MASS_MEV;

	// Create parameter struct
	containers::parameters params{V0, kappa, lambda, r0, a, Rls, als, N, Z, l, k, xmin, xmax,
		xmatch, m, tensorV, kappa_so, scenario, isospin};
	
	// Run dirac solver
	std::pair<double, double> result = solveDirac(params, a0, B);  
	double B_result = result.first, a0_result = result.second;
	std::cout << "Converged values are B: " << B_result << ", and a0: " << a0_result;
	return 0;
}
