/*
Project:        Shell evolution of the dirac equation
                
Authors:        Alexander Kiessling
                (2022-2023)

Description:    Main interface file of program, initializes data
                and calls solver routines.
*/

#include "solver.h"
#include "util.h"

int main()
{
	// Load data from files, get user input
	containers::parameters params = setup();
	
	// Run dirac solver
	containers::solutionsBa0 result = solveDirac(params, params.a0, params.B);  
	
	// print and save wavefunction to txt file for plotting
	double B_result = result.B, a0_result = result.a0;
	std::cout << "Converged values are B: " << B_result << ", and a0: " << a0_result << "\n";
	saveWF(result.rvals, result.fvals, result.gvals);
	
	
	return 0;
}
