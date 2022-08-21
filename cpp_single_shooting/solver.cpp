#include "solver.h"

std::pair <double, double> Tensor_fit(double x, double F, double G, const containers::parameters& params)
{
	double Col = 0;

	if (params.isospin == 1)
	{
		if (x > params.sigmaR)
		{
			Col = 0.0072923 * params.Z / x;
		}
		else
		{
			Col = 0.0072923 * params.Z * (3 * pow(params.sigmaR, 2) - pow(x, 2)) / (2 * pow(params.sigmaR, 3));
		}
	}
	double sigma = params.sigmaV0 / (1 + exp((x - params.sigmaR) / params.sigmaa)) + Col;
	double delta = params.dV0 / (1 + exp((x - params.dR) / params.da)) + Col;
	double U = params.tensorV / (1 + exp((x - params.sigmaR) / params.sigmaa));

	double dfgW_1 = (-params.B + sigma) * G + (params.k / x - U) * F;
	double dfgW_2 = (2 * params.m + params.B - delta) * F + (U - params.k / x) * G;
	return { dfgW_1, dfgW_2 };
}

std::pair <double, double> int_n_Tensor(double xstart, double xend, double iniF, double iniG,
	const containers::parameters& params)
{
	double h = 0.001;
	double step = (xend - xstart) * h;
	double F = iniF;
	double G = iniG;
	double step_2 = step / 2;
	double step_6 = step / 6;
	int end_condition = int(1.0 / h);

	double x = 0;

	std::pair <double, double> v1, v2, v3, v4;

	for (int i = 0; i++; i < end_condition)
	{
		x = i * step + xstart;
		v1 = Tensor_fit(x, F, G, params);
		v2 = Tensor_fit(x + step_2, F + v1.first * step_2, G + v1.second * step_2, params);
		v3 = Tensor_fit(x + step_2, F + v2.first * step_2, G + v2.second * step_2, params);
		v4 = Tensor_fit(x + step, F + v3.first * step, G + v3.second * step, params);
		F += (v1.first + 2 * v2.first + 2 * v3.first + v4.first) * step_6;
		G += (v1.second + 2 * v2.second + 2 * v3.second + v4.second) * step_6;
	}

	return { F, G };
}

std::tuple<int, double, double> spectr(std::string state)
{
	std::istringstream ss(state);

	std::string word;
	std::vector<std::string> sentence;
	while (ss >> word)
	{
		sentence.push_back(word);
	}

	std::map<char, int> level_map = {
	{ 's', 0 },
	{ 'p', 1 },
	{ 'd', 2 },
	{ 'f', 3 },
	{ 'g', 4 },
	{ 'h', 5 },
	{ 'i', 6 },
	{ 'j', 7 },
	{ 'k', 8 }
	};

	int l = level_map[sentence[1][1]];
	double j = (double(sentence[1][2]) - 48) / 2;
	double k;

	if (j == l + 1 / 2)
	{
		k = -(l + 1);
	}
	else
	{
		k = l;
	}
	return { l, j, k };
}

std::tuple<double, double, double, double> BC(double xmin, double xmax, const containers::parameters& params)
{
	double Foutbc, Goutbc, miu, Finbc, Ginbc;
	if (params.k < 0)
	{
		Foutbc = -params.a0 * pow(xmin, (params.l + 2)) * (-params.B + params.sigmaV0 / (1 + exp(-params.sigmaR / params.sigmaa))) / (params.l + 2 - params.k);
		Goutbc = params.a0 * pow(xmin, (params.l + 1));
	}
	else
	{
		Foutbc = params.a0 * pow(xmin, params.l);
		Goutbc = params.a0 * pow(xmin, (params.l + 1)) * (2 * params.m + params.B - params.dV0 / (1 + exp(-params.dR / params.da))) / (params.l + params.k + 1);

	}
	try
	{
		miu = sqrt(-2 * params.m * params.B - pow(params.B, 2));
		Finbc = -sqrt(-params.B / (2 * params.m + params.B)) * exp(-miu * xmax);
	}
	catch (...)
	{
		miu = 1e-10;
		Finbc = -sqrt(-params.B / (1e-10)) * exp(-miu * xmax);
	}
	Ginbc = exp(-miu * xmax);

	return { Foutbc, Goutbc, Finbc, Ginbc };
}

std::tuple<double, double, double, double> BC_pos(double xmin, double xmax, const containers::parameters& params)
{
	double Foutbc, Goutbc, Finbc, Ginbc;
	if (params.k < 0)
	{
		Foutbc = -params.a0 * pow(xmin, (params.l + 2)) * (-params.B + params.sigmaV0 / (1 + exp(-params.sigmaR / params.sigmaa))) / (params.l + 2 - params.k);
		Goutbc = params.a0 * pow(xmin, (params.l + 1));
		Finbc = (std::sph_bessel(params.l, xmax) + std::sph_neumann(params.l, xmax));
		Ginbc = sqrt(pow(params.B, 2) + 2 * params.B * params.m) / (params.B + 2 * params.m) * (std::sph_bessel(params.l + 1, xmax) + std::sph_neumann(params.l + 1, xmax));
	}
	else
	{
		Foutbc = params.a0 * pow(xmin, params.l);
		Goutbc = params.a0 * pow(xmin, (params.l + 1)) * (2 * params.m + params.B - params.dV0 / (1 + exp(-params.dR / params.da))) / (params.l + params.k + 1);
		Finbc = (std::sph_bessel(params.l, xmax) + std::sph_neumann(params.l, xmax));
		Ginbc = sqrt(pow(params.B, 2) + 2 * params.B * params.m) / (params.B + 2 * params.m) * (std::sph_bessel(params.l - 1, xmax) + std::sph_neumann(params.l - 1, xmax));
	}
	double norm = sqrt(pow(std::sph_bessel(params.l, xmax), 2) + pow(std::sph_neumann(params.l, xmax), 2));

	return { Foutbc, Goutbc, Finbc / norm, Ginbc / norm };
}

/* RETURN TO GAP FUNCTION LATER
double Gap(x, y) :
	return np.array(x - y).T

*/

bool test_converge(double B, double Exp)
{
	if (Exp < 0)
	{
		if (B > Exp - 20 && B < 0)
		{
			return true;
		}
	}
	else if (Exp > 0)
	{
		if (B < Exp + 20 && B > 0)
		{
			return true;
		}
	}
	return false;
}

void solve_dirac()
{

}

int main()
{
	int element = 0;
	int scenario = 1;

	// File pointer
	fstream fin;

	// Open an existing file
	fin.open("reportcard.csv", ios::in);


	std::cout << "hi";
	return 0;
}
