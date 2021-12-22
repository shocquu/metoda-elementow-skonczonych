#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>

#include "../utils/gauss.h"
#include "../utils/matrix.h"
#include "../utils/helpers.h"
#include "element4_2d.h"
#include "element.h"
#include "node.h"

static double alpha = 300;
static double simTime = 100;
static double simStepTime = 50;
static double conductivity = 25;
static double initialTemp = 100;
static double ambientTemp = 1200;
static double specificHeat = 700;
static double density = 7800;

struct Grid {
	Element* elements;
	Node* nodes;
	double width, height;
	int nW, nH, nN, nE, p;

	double** aggrH, ** aggrC;
	double* aggrP;

	Grid();
	Grid(double height, double width, int nH, int nW);
	~Grid();
		
	void launch(Element4_2D &el4, std::string path, bool showDetails = false);
	void launch(Element4_2D &el4, double simTime, double simStepTime, double conductivity, double alpha, double density, double specificHeat, double initialTemp, double ambientTemp, bool showDetails = false);
	void readFromFile(std::string path, double& simulationTime, double& simulationStepTime, double& conductivity, double& alpha, double& density, double& specificHeat, double& initialTemp, double& ambientTemp);
	void collectData(Element4_2D el4);
	void aggregate(Element &currEl);
	void calcHbc(Element4_2D &el4, Element &currEl, double alpha = alpha, double ambientTemp = ambientTemp);
	void calcH(Element4_2D &el4, Element &currEl, double k = conductivity);
	void calcC(Element4_2D &el4, Element &currEl, double c = specificHeat, double ro = density);
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
};