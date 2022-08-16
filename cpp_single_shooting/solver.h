#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <tuple>
#include <iostream>

#define NEUTRON_MASS_MEV 939.56542052
#define PROTON_MASS_MEV 938.27208816

namespace containers
{
	typedef struct
	{
		double a0;
		double B;
		double sigmaV0;
		double sigmaR;
		double sigmaa;
		double dV0;
		double dR;
		double da;
		double k;
		double m;
		double Z;
		double tensorV;
		int isospin;
		int l;

	} parameters;
}