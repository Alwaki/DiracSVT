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
        double V0 = 0;
        double kappa = 0;
        double lambda = 0;
        double r0 = 0;
        double a = 0;
        double Rls = 0;
        double als = 0;
        double isospin = 0;
        double N = 0;
        double Z = 0;
        double l = 0;
        double k = 0;
        double xmin = 0;
        double xmax = 0;
        double xmatch = 0;
        double m = 0;
        int scenario = 0;
        double tensorV = 0;

        parameters(double v0, double v1, double v2, double v3, 
            double v4, double v5, double v6, double v7, double v8,
            double v9, double v10, double v11, double v12, double v13, 
            double v14, double v15, int v16, double v17)
        {
            double V0 = v0;
            double kappa = v1;
            double lambda = v2;
            double r0 = v3;
            double a = v4;
            double Rls = v5;
            double als = v6;
            double isospin = v7;
            double N = v8;
            double Z = v9;
            double l = v10;
            double k = v11;
            double xmin = v12;
            double xmax = v13;
            double xmatch = v14;
            double m = v15;
            int scenario = v16;
            double tensorV = v17;
        }

	} parameters;
}