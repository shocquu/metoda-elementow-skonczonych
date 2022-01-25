#pragma once

#include <math.h>
#include <vector>
#include <array>

#include "../utils/gauss.h"
#include "element.h"

struct Side { double* pc1, * pc2; };

struct Element4_2D {
	double dNdX[4][4] = { 0 }, dNdY[4][4] = { 0 }, Npc1[4][4] = { 0 }, Npc2[4][4] = { 0 };
	std::vector<std::array<double, 4>> nMatrix, dNdKsiMatrix, dNdEtaMatrix;
	std::vector<double> ksi, eta, w;
	double detJ;
	int p, n;

	Element4_2D(int integralPoints = 2);
	
	std::array<double, 4> nKsiEta(double ksi, double eta);
	std::array<double, 4> dNdEta(double ksi);
	std::array<double, 4> dNdKsi(double eta);
};