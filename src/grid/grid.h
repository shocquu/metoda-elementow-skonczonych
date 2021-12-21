#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include "node.h"
#include "element.h"

struct Grid {
	Element* elements;
	Node* nodes;
	double width, height;
	int nW, nH, nN, nE;

	Grid();
	Grid(double height, double width, int nH, int nW);
	~Grid();

	void readFromFile(std::string path, double& simulationTime, double& simulationStepTime, double& conductivity, double& alpha, double& density, double& specificHeat, double& initialTemp, double& ambientTemp);
	void setNodesCoords();
	void printNodesCoords();
	void setElementsNodes();
	void printElementsNodes();
	void setNodesBC();
	void printBCNodes();
};