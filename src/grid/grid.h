#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <tuple>

#include "../../include/matplotlib-cpp/matplotlibcpp.h"
#include "../../include/meshgrid.hpp"
#include "../utils/gauss.h"
#include "../utils/matrix.h"
#include "../utils/helpers.h"
#include "element4_2d.h"
#include "element.h"
#include "node.h"

struct Grid {
	Element* elements;
	Node* nodes;
	double width, height;
	int nW, nH, nN, nE, p;
	double** aggrH, ** aggrC;
	double* aggrP;

	Grid();
	Grid(Element4_2D el4, std::string path);
	Grid(double height, double width, int nH, int nW);
	~Grid();
		
	void plotHeatMap(Element4_2D el4);
	void calcTemperature(Element4_2D el4, bool showMinMax = false, bool saveToFile = false);
	void launch(Element4_2D &el4, std::string path, bool showDetails = false);
	void launch(Element4_2D &el4, double simTime, double simStepTime, double conductivity, double alpha, double density, double specificHeat, double initialTemp, double ambientTemp, bool showDetails = false);
	void readFromFile(std::string path);
	void collectData(Element4_2D el4);
	void aggregate(Element &currEl);
	void calcHbc(Element4_2D &el4, Element &currEl, double alpha, double ambientTemp);
	void calcH(Element4_2D &el4, Element &currEl, double k);
	void calcC(Element4_2D &el4, Element &currEl, double c, double ro);
	void calcJ(Element4_2D &el4, Element &currEl);
	void jacobian(Element4_2D &el4, Element &currEl);
	double distance(double x1, double y1, double x2, double y2);
	void initMatrices();
	void setNodesCoords();
	void printNodesCoords();
	void setElementsNodes();
	void printElementsNodes();
	void setNodesBC();
	void printBCNodes();
	void printSimulationData();

private:
	static double simTime, simStepTime, conductivity, alpha, density, specificHeat, initialTemp, ambientTemp;
};

