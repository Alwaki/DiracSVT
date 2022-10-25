#include "solver.h"
#include <Eigen/Core>

/* RETURN TO GAP FUNCTION LATER
double Gap(x, y) :
	return np.array(x - y).T

*/
/*
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
*/

std::pair <double, double> tensorFit(double x, double F, double G, const containers::parameters& params, 
										double B, double sigmaV0, double sigmaR, double sigmaa, double dV0, 
										double dR, double da)
{
	double Col = 0;

	if (params.isospin == 1)
	{
		if (x > sigmaR)
		{
			Col = 0.0072923 * params.Z / x;
		}
		else
		{
			Col = 0.0072923 * params.Z * (3 * pow(sigmaR, 2) - pow(x, 2)) / (2 * pow(sigmaR, 3));
		}
	}
	
	double sigma = sigmaV0 / (1 + exp((x - sigmaR) / sigmaa)) + Col;
	double delta = dV0 / (1 + exp((x - dR) / da)) + Col;
	double U = params.tensorV / (1 + exp((x - sigmaR) / sigmaa));
	double dfgW_1 = (-B + sigma) * G + (params.k / x - U) * F;
	double dfgW_2 = (2 * params.m + B - delta) * F + (U - params.k / x) * G;
	return { dfgW_1, dfgW_2 };
}

std::pair <double, double> integrateTensor(double iniF, double iniG, const containers::parameters& params, 
											double B, double a0, double sigmaV0, double sigmaR, double sigmaa, 
											double dV0, double dR, double da, double xend, double xstart)
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

	for (int i = 0; i < end_condition; i++)
	{
		x = i * step + xstart;
		v1 = tensorFit(x, F, G, params, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		v2 = tensorFit(x + step_2, F + v1.first * step_2, G + v1.second * step_2, 
			 params, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		v3 = tensorFit(x + step_2, F + v2.first * step_2, G + v2.second * step_2, 
			 params, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		v4 = tensorFit(x + step, F + v3.first * step, G + v3.second * step, 
			params, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		F += (v1.first + 2 * v2.first + 2 * v3.first + v4.first) * step_6;
		G += (v1.second + 2 * v2.second + 2 * v3.second + v4.second) * step_6;
	}

	return { F, G };
}

std::tuple<double, double, double, double> BC(const containers::parameters& params, double B, double a0, double sigmaV0,
											double sigmaR, double sigmaa, double dV0, double dR, double da)
{
	double Foutbc, Goutbc, miu, Finbc, Ginbc;
	if (params.k < 0)
	{
		Foutbc = -a0 * pow(params.xmin, (params.l + 2)) * (-B + sigmaV0 / (1 + exp(-sigmaR / sigmaa))) / (params.l + 2 - params.k);
		Goutbc = a0 * pow(params.xmin, (params.l + 1));
	}
	else
	{
		Foutbc = a0 * pow(params.xmin, params.l);
		Goutbc = a0 * pow(params.xmin, (params.l + 1)) * (2 * params.m + B - dV0 / (1 + exp(-dR / da))) / (params.l + params.k + 1);

	}
	try
	{
		miu = sqrt(-2 * params.m * B - pow(B, 2));
		Finbc = -sqrt(-B / (2 * params.m + B)) * exp(-miu * params.xmax);
	}
	catch (...)
	{
		miu = 1e-10;
		Finbc = -sqrt(-B / (1e-10)) * exp(-miu * params.xmax);
	}
	Ginbc = exp(-miu * params.xmax);

	return { Foutbc, Goutbc, Finbc, Ginbc };
}

std::tuple<double, double, double, double> BC_pos(const containers::parameters& params, double B, double a0, double sigmaV0,
											double sigmaR, double sigmaa, double dV0, double dR, double da)
{
	double Foutbc, Goutbc, Finbc, Ginbc;
	if (params.k < 0)
	{
		Foutbc = -a0 * pow(params.xmin, (params.l + 2)) * (-B + sigmaV0 / (1 + exp(-sigmaR / sigmaa))) / (params.l + 2 - params.k);
		Goutbc = a0 * pow(params.xmin, (params.l + 1));
		Finbc = (std::sph_bessel(params.l, params.xmax) + std::sph_neumann(params.l, params.xmax));
		Ginbc = sqrt(pow(B, 2) + 2 * B * params.m) / (B + 2 * params.m) * (std::sph_bessel(params.l + 1, params.xmax) + std::sph_neumann(params.l + 1, params.xmax));
	}
	else
	{
		Foutbc = a0 * pow(params.xmin, params.l);
		Goutbc = a0 * pow(params.xmin, (params.l + 1)) * (2 * params.m + B - dV0 / (1 + exp(-dR / da))) / (params.l + params.k + 1);
		Finbc = (std::sph_bessel(params.l, params.xmax) + std::sph_neumann(params.l, params.xmax));
		Ginbc = sqrt(pow(B, 2) + 2 * B * params.m) / (B + 2 * params.m) * (std::sph_bessel(params.l - 1, params.xmax) + std::sph_neumann(params.l - 1, params.xmax));
	}
	double norm = sqrt(pow(std::sph_bessel(params.l, params.xmax), 2) + pow(std::sph_neumann(params.l, params.xmax), 2));

	return { Foutbc, Goutbc, Finbc / norm, Ginbc / norm };
}

std::pair <double, double> solveDirac(const containers::parameters& params, double a0_in, double B0)
{
	// Setup potential
	double A = params.N + params.Z;
	double sigmaV0 = 0;
	double dV0 = 0;
	if (params.isospin == 1)
	{
		sigmaV0 = params.V0*(1 + params.kappa*(params.N-params.Z)/A);
		if (params.scenario == 1)
		{
			dV0 = -params.lambda*sigmaV0;
		}
		else if (params.scenario==2)
		{
			dV0 = -params.lambda*params.V0*(1 - params.kappa*(params.N-params.Z)/A);
		}
		else if (params.scenario==3)
		{
			dV0 = -params.lambda*params.V0*(1 - params.kappa_so*(params.N-params.Z)/A);
		}
	}
	else if (params.isospin == -1)
	{
		sigmaV0 = params.V0*(1 - params.kappa*(params.N-params.Z)/A);
		if (params.scenario==1)
		{
			dV0 = -params.lambda * sigmaV0;
		}
		else if (params.scenario==2)
		{
			dV0 = -params.lambda*params.V0*(1 + params.kappa*(params.N-params.Z)/A);
		}
		else if (params.scenario==3)
		{
			dV0 = -params.lambda*params.V0*(1 + params.kappa_so*(params.N-params.Z)/A);
		}
		
		
	}
    double sigmaR = params.r0*pow(A,1.0/3.0);
    double dR = params.Rls*pow(A,1.0/3.0);
    double sigmaa = params.a;
    double da = params.a;
    double error = 100;
	double B = B0;
	double a0 = a0_in;

	// Iterate solvers
	double h = 0.0001;
	int iterations = 0;
	while (error > 0.0001)
	{
		double Foutbc, Goutbc, Finbc, Ginbc;
		if(B <0)
		{
			std::tie(Foutbc, Goutbc, Finbc, Ginbc) = BC(params, B, a0, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		}
		else
		{
			std::tie(Foutbc, Goutbc, Finbc, Ginbc) = BC_pos(params, B, a0, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		}
		std::pair <double, double> inFG = integrateTensor(Finbc, Ginbc, params, B, a0, sigmaV0, 
										sigmaR, sigmaa, dV0, dR, da, params.xmatch, params.xmax);
        std::pair <double, double> outFG = integrateTensor(Foutbc, Goutbc, params, B, a0, sigmaV0, 
										sigmaR, sigmaa, dV0, dR, da, params.xmatch, params.xmin);
		double B1 = B + B*h;
		
		if(B1 < 0)
		{
			std::tie(Foutbc, Goutbc, Finbc, Ginbc) = BC(params, B1, a0, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		}
		else
		{
			std::tie(Foutbc, Goutbc, Finbc, Ginbc) = BC_pos(params, B1, a0, sigmaV0, sigmaR, sigmaa, dV0, dR, da);
		}
		std::pair <double, double> dBinFG = integrateTensor(Finbc, Ginbc, params, B1, a0, sigmaV0, 
										sigmaR, sigmaa, dV0, dR, da, params.xmatch, params.xmax);
        std::pair <double, double> dBoutFG = integrateTensor(Foutbc, Goutbc, params, B1, a0, sigmaV0, 
										sigmaR, sigmaa, dV0, dR, da, params.xmatch, params.xmin);
										
		double dGFdB_1    = ((dBoutFG.first - dBinFG.first) - (outFG.first - inFG.first))/(B*h);
		double dGFdB_2    = ((dBoutFG.second - dBinFG.second) - (outFG.second - inFG.second))/(B*h);
		double da0outGF_1 = outFG.first * (1.0 + h);
		double da0outGF_2 = outFG.second * (1.0 + h);
		double dGFda0_1   = ((da0outGF_1 - inFG.first) - (outFG.first - inFG.first)) / (a0*h);
		double dGFda0_2   = ((da0outGF_2 - inFG.second) - (outFG.second - inFG.second)) / (a0*h);
		double dOutIn_1   = (outFG.first - inFG.first);
		double dOutIn_2   = (outFG.second - inFG.second);

		Eigen::Matrix2f M;
		M <<  dGFdB_1, dGFdB_2, dGFda0_1, dGFda0_2;
		Eigen::Vector2f Old, Cold, diff, New;
		Old << B, a0;
		diff << dOutIn_1, dOutIn_2;
		Cold = (M*Old)-diff;
		New = M.completeOrthogonalDecomposition().solve(Cold);
		B = New(1);
		a0 = New(2);
		error = (New-Old).norm();

		iterations++;
	}

return {B, a0};
}

int read_user_input(std::string type)
{
	int temp;
	std::cout << "Enter " << type << ": ";
	std::cin >> temp;
	if (type == "scenario")
	{
		while (1)
		{
			if (std::cin.fail() || (temp < 1) || (temp > 3))
			{
				std::cin.clear();
				std::cin.ignore(1000, '\n');
				std::cout << "Invalid input, please try again. " << std::endl;
				std::cin >> temp;
			}
			else
			{
				break;
			}
		}
	}
	else
	{
		while (1)
		{
			if (std::cin.fail() || ( temp < 1) || (temp > 89) )
			{
				std::cin.clear();
				std::cin.ignore(1000, '\n');
				std::cout << "Invalid input, please try again. " << std::endl;
				std::cin >> temp;
			}
			else
			{
				break;
			}
		}
	}
	return temp;
}

std::vector<std::vector<double>> read_csv(std::string filename) {
	// Reads a CSV file into a vector of <string, vector<int>> pairs where
	// each pair represents <column name, column values>

	// Create an input filestream
	std::ifstream myFile(filename);

	// Make sure the file is open
	if (!myFile.is_open()) throw std::runtime_error("Could not open file");

	std::vector<std::vector<double>> content;
	std::vector<double> row;
	std::string line, word;

	while (getline(myFile, line))
	{
		row.clear();

		std::stringstream str(line);

		while (getline(str, word, ','))
		{
			try
			{
				row.push_back(std::stod(word));
			}
			catch (...)
			{

			}
		}
		content.push_back(row);
	}
	/*
	for (int i = 0; i < content.size(); i++)
	{
		for (int j = 0; j < content[i].size(); j++)
		{
			std::cout << content[i][j] << " ";
		}
		std::cout << "\n";
	} */

	return content;
}

int main()
{
	// Allow user to select scenario and state
	std::cout << "The program requires a scenario and state. These \n";
	std::cout << "are both integers, and the scenario ranges from 1 to 3 \n";
	std::cout << "while the state ranges from 1 to 89. \n";
	std::cout << "\n";
	int scenario = 1;//read_user_input("scenario"); //TODO: Remove these comments. However, for debug keep.
	int state = 1;//read_user_input("state");

	// Load data from files
	auto data = read_csv("data.csv");
	auto parameters = read_csv("parameters.csv");

	// Set parameters dependent on scenario
	double V0 = parameters[scenario][1];
	double kappa = parameters[scenario][2];
	double lambda = parameters[scenario][3];
	double r0 = parameters[scenario][4];
	double a = parameters[scenario][5];
	double Rls = parameters[scenario][6];
	double als = parameters[scenario][7];

	// Set parameters dependent on state
	int isospin = data[state][2];
	double N = data[state][3];
	double Z = data[state][4];
	double l = data[state][5];
	double k = data[state][6];
	double B = data[state][7];
	double xmin = data[state][8];
	double xmax = data[state][9];
	double xmatch = data[state][10];
	double a0 = data[state][11];

	// Set other parameters
	double tensorV = 0;
	double kappa_so = 0;
	double m = 0;
	if (isospin == -1)
		m = NEUTRON_MASS_MEV;
	else
		m = PROTON_MASS_MEV;

	// Create parameter struct
	containers::parameters params{V0, kappa, lambda, r0, a, Rls, als, N, Z, l, k, xmin, xmax,
		xmatch, m, tensorV, kappa_so, scenario, isospin};
	
	// Run dirac solver
	std::pair<double, double> result = solveDirac(params, a0, B);  
	double B_result = result.first, a0_result = result.second;
	std::cout << "Converged values are B: " << B_result << ", and a0: " << a0_result;
	return 0;
}
