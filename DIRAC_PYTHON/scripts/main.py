#!/usr/bin/env python3

"""
Title:        		DiracSVT
                
Authors:        	Alexander Kiessling, Daniel Karlsson, 
			Yuxin Zhao, Chong Qi

Version:		1.1 (11/2023)	

Project Description:    Numerical solution of the Dirac equation with scalar,
			vector and tensor potentials

File Description:	Main interface file of program, parses user arguments
                        and calls solver routines.
"""

from util import *
from shooting_method import *

def main():
    parser = argparse.ArgumentParser(description="Dirac solver")
    parser.add_argument(
            "--state",
            default=1,
            type=int,
            help="Integer value from 1 to 89 for specific state, 0 runs all."
    )
    parser.add_argument(
            "--scenario",
            default=1,
            type=int,
            choices=[1,2,3],
            help="Integer value from 1 to 3, specifies parameters and potential"
    )

    args = parser.parse_args()
    
    if args.state == 0:
        for state in range(1,89):
                # Run solver routine
                B, a0, rvals, FGvals = run_1Dsolve(state, args.scenario)
                print("B: ", B, " a0: ", a0)  
                
                
    else:
        # Run solver routine
        B, a0, rvals, FGvals = run_1Dsolve(args.state, args.scenario)
        
        # If plotting enabled, plot wavefunction. Print results regardless.
        print("B: ", B, " a0: ", a0)
        plotWF(rvals, FGvals)


    return 0

if __name__ == "__main__":
    sys.exit(main())
