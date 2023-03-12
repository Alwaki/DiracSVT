/*
Title:        			DiracSVT
                
Authors:        		Alexander Kiessling, Daniel Karlsson, 
						Yuxin Zhao, Chong Qi

Version:				1.0 (03/2023)	

Project Description:    Numerical solution of the Dirac equation with scalar,
						vector and tensor potentials

File Description:		Main interface file of program, initializes data
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
		for(int iter = 1; iter < 90; iter++)
		{
			// Load data from files
			containers::parameters params = setup(selection.first, iter);
			
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
		containers::parameters params = setup(selection.first, selection.second);
		
		// Run dirac solver
		containers::solutionsBa0 result = solveDirac(params, params.a0, params.B);  
		
		// print and save wavefunction to txt file for plotting
		double B_result = result.B, a0_result = result.a0;
		std::cout << "Converged values are B: " << B_result << ", and a0: " << a0_result << "\n";
		saveWF(result.rvals, result.fvals, result.gvals);
	}
	return 0;
}
