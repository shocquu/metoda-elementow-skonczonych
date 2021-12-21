#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>

#include "../utils/matrix.h"
#include "element4_2d.h"
#include "element.h"
#include "node.h"

struct Grid {
	Element* elements;
	Node* nodes;
	double width, height;
	int nW, nH, nN, nE;

	double* globalP;
	double** globalH, ** globalC;

	Grid();
	Grid(double height, double width, int nH, int nW);
	~Grid();

	void initMatrices();
	void readFromFile(std::string path, double& simulationTime, double& simulationStepTime, double& conductivity, double& alpha, double& density, double& specificHeat, double& initialTemp, double& ambientTemp);
	void aggregate(Element currEl);
	void calcHbc(Element4_2D &el4, Element &currEl, int elIndex, double alpha, double ambientTemp);
	void calcH(Element4_2D &el4, Element &currEl, double k);
	void calcC(Element4_2D &el4, Element &currEl, double c, double ro);
	void calcJ(Element4_2D &el4, double J[2][2], int pc, double x, double y);
	void jacobian(Element4_2D &el4, Element &currEl);
	double distance(double x1, double y1, double x2, double y2);
	void setNodesCoords();
	void printNodesCoords();
	void setElementsNodes();
	void printElementsNodes();
	void setNodesBC();
	void printBCNodes();
};