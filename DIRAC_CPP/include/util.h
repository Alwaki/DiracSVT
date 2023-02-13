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

#endif