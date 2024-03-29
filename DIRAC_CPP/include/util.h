/*
Title:        			DiracSVT
                
Authors:        		Alexander Kiessling, Daniel Karlsson, 
						Yuxin Zhao, Chong Qi

Version:				1.1 (11/2023)	

Project Description:    Numerical solution of the Dirac equation with scalar,
						vector and tensor potentials

File Description:		Header file for util.cpp. Also serves as 
                        function documentation.
*/


#ifndef UTIL
#define UTIL

#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <tuple>
#include <iostream>
#include <fstream>
#include <limits>
#include <utility>

#include "solver.h"

/*
* brief: Takes in user input in the form of scenario and state selection
*/
int read_user_input(std::string type);

/*
* brief: Loads data from csv files to get parameters
*/
std::vector<std::vector<double>> read_csv(std::string filename);

/*
* brief: Saves solution of wavefunction to .txt file
*/
void saveWF(std::vector<double> rvals, std::vector<double> fvals, std::vector<double> gvals);

/*
* brief: Convenience function to load data and collect it in struct 
*/
containers::parameters setup(int scenario, int state);

/*
* brief: Convenience function to get user input
*/
std::pair<int,int> user_selection();

#endif
