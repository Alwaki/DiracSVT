/*
Project:        Shell evolution of the dirac equation
                
Authors:        Alexander Kiessling
                (2022-2023)

Description:    Contains all utility function definitioons such as reading 
				and loading input and data. See header file for more 
				information.
*/


#include "util.h"

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
			if (std::cin.fail() || ( temp < 0) || (temp > 89) )
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

	return content;
}

void saveWF(std::vector<double> rvals, std::vector<double> fvals, std::vector<double> gvals)
{
	std::ofstream fw("wavefunction.txt", std::ofstream::out);
	if (fw.is_open())
    {
      //store array contents to text file
      for (int i = 0; i < rvals.size(); i++) 
	  {
        fw << rvals[i] << ", " << fvals[i] << ", " << gvals[i] << "\n";
      }
      fw.close();
    }
}

std::pair<int,int> user_selection()
{
	// Allow user to select scenario and state
	std::cout << "The program requires a scenario and state. These \n";
	std::cout << "are both integers, and the scenario ranges from 1 to 3 \n";
	std::cout << "while the state ranges from 1 to 89. If all states \n";
	std::cout << "are to be run, select state 0. \n";
	std::cout << "\n";
	int scenario = read_user_input("scenario"); 
	int state = read_user_input("state");
	return {scenario, state};
}

containers::parameters setup(int scenario, int state)
{
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
		xmatch, m, tensorV, kappa_so, B, a0, scenario, isospin};
    
    return params;
}