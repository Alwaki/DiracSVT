#!/usr/bin/env python3

"""
Project:        Shell evolution of the dirac equation
                
Authors:        Daniel Karlsson & Alexander Kiessling
                (2021-2023)

Description:    Main interface file of program, parses user arguments
                and calls solver routines.
"""

from util import *
from shooting_method import *

def main():
    parser = argparse.ArgumentParser(description="Dirac solver")
    parser.add_argument(
            "--state",
            default='16O 1p1/2',
            help="The atom and state to investigate, ex: '16O 1p1/2'"
    )
    parser.add_argument(
            "--particle",
            default=-1,
            type=int,
            choices=[1,-1],
            help="Sets particle type, -1 for neutron, 1 for proton."
    )
    parser.add_argument(
            "--scenario",
            default=1,
            type=int,
            choices=[1,2,3],
            help="Integer value from 1 to 3, specifies parameters and potential"
    )
    parser.add_argument(
            "--plot",
            default=False,
            help="Boolean for if plotting wavefunction should be enabled"
    )

    args = parser.parse_args()
    
    # Create state
    state = [args.state, args.particle]
    
    # Run solver routine
    B, a0, rvals, FGvals = run_1Dsolve(state, Scenario=args.scenario)
    
    # If plotting enabled, plot wavefunction. Print results regardless.
    print("B: ", B, " a0: ", a0)
    if args.plot==True:
        print(args.plot)
        plotWF(rvals, FGvals, state)

    return 0

if __name__ == "__main__":
    sys.exit(main())
