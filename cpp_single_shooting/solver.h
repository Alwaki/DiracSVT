#pragma once

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
        double V0, kappa, lambda, r0, a, Rls, als, isospin, N, Z, l, k, xmin, xmax, xmatch,
                m, tensorV, kappa_so;
        int scenario;

	} parameters;
}