#ifndef UTIL
#define UTIL

//#include "matplotlibcpp.h"
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

//namespace plt = matplotlibcpp;


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

bool plotWF(std::vector<double> rvals, std::vector<double> fvals, std::vector<double> gvals)
{
	/*
	plt::figure();
	plt::plot(rvals, fvals, "b");
	plt::plot(rvals, gvals, "r");
	plt::show();
	*/
	/*
	RGBABitmapImageReference *imageReference = CreateRGBABitmapImageReference();
	StringReference *errorMessage = CreateStringReferenceLengthValue(0, L' ');

	DrawScatterPlot(imageReference, 600, 400, &rvals, &fvals, errorMessage);
	DrawScatterPlot(imageReference, 600, 400, &rvals, &gvals, errorMessage);
	vector<double> *pngdata = ConvertToPNG(imageReference->image);
	WriteToFile(pngdata, "plotWF.png");
	DeleteImage(imageReference->image);
	FreeAllocations();
	*/
}

#endif