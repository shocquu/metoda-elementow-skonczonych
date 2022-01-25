#pragma once

#include <math.h>
#include <vector>
#include <array>

#include "../utils/gauss.h"
#include "element.h"

struct Element4_2D {
	// double **nMatrix, **dNdKsiMatrix, **dNdEtaMatrix;
	double dNdX[4][4] = { 0 }, dNdY[4][4] = { 0 };
	std::vector<double> ksi, eta;
	std::vector<std::array<double, 4>> nMatrix, dNdKsiMatrix, dNdEtaMatrix;
	double detJ;
	double *w;
	int p;

	Element4_2D(int integralPoints = 2);
	
	std::array<double, 4> nKsiEta(double ksi, double eta);
	std::array<double, 4> dNdEta(double ksi);
	std::array<double, 4> dNdKsi(double eta);
};