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
	// get user input
	std::pair selection = user_selection();

	if(selection.second == 0)
	{
		for(int i = 1; i < 90; i++)
		{
			// Load data from files
			containers::parameters params = setup();
			
			// Run dirac solver
			containers::solutionsBa0 result = solveDirac(params, params.a0, params.B);  
			
			// print converged values
			double B_result = result.B, a0_result = result.a0;
			std::cout << "Converged values are B: " << B_result << ", and a0: " << a0_result << "\n";
		}
	}
	else
	{
		// Load data from files
		containers::parameters params = setup();
		
		// Run dirac solver
		containers::solutionsBa0 result = solveDirac(params, params.a0, params.B);  
		
		// print and save wavefunction to txt file for plotting
		double B_result = result.B, a0_result = result.a0;
		std::cout << "Converged values are B: " << B_result << ", and a0: " << a0_result << "\n";
		saveWF(result.rvals, result.fvals, result.gvals);
	}
	return 0;
}
