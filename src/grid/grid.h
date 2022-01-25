#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <tuple>
#include <array>
#include <vector>

#include "../../include/matplotlib-cpp/matplotlibcpp.h"
#include "../../include/meshgrid.hpp"
#include "../utils/gauss.h"
#include "../utils/matrix.h"
#include "../utils/helpers.h"
#include "element4_2d.h"
#include "element.h"
#include "node.h"

struct SimulationData {
	double alpha = 300;
	double simTime = 60;
	double simStepTime = 1;
	double conductivity = 25;
	double density = 7800;
	double specificHeat = 700;
	double initialTemp = 100;
	double ambientTemp = 1200;
};

struct Grid {
	SimulationData data;
	Element4_2D el4;
	Element* elements;
	Node* nodes;
	double width, height;
	int nW, nH, nN, nE, p;
	double** aggrH, ** aggrC;
	double* aggrP;

	Grid();
	Grid(std::string path);
	Grid(double height, double width, int nH, int nW);
	~Grid();
		
	Grid& heatMap();
	Grid& start(bool saveToFile = false);
	void plotHeatMap();
	void launch(std::string path, bool saveToFile = false);
	void launch(SimulationData data, bool saveToFile = false);
	void readFromFile(std::string path);
	void aggregate(Element &currEl);
	void calcHbc(Element &currEl);
	void calcH(Element &currEl);
	void calcC(Element &currEl);
	void calcJ(Element &currEl, int pc);
	void jacobian(Element &currEl);
	void calcTemperature(bool showMinMax = false, bool saveToFile = false);
	double distance(double x1, double y1, double x2, double y2);
	void initMatrices();
	void setNodesCoords();
	void printNodesCoords();
	void setElementsNodes();
	void printElementsNodes();
	void setNodesBC();
	void printBCNodes();
	void printSimulationData();
};

