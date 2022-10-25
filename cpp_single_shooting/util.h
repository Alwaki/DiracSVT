#pragma once
#include "solver.h"


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