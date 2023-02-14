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

int read_user_input(std::string type);

std::vector<std::vector<double>> read_csv(std::string filename);

void saveWF(std::vector<double> rvals, std::vector<double> fvals, std::vector<double> gvals);

containers::parameters setup();

#endif