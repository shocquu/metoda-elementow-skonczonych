#pragma once

#include <math.h>

#include "../utils/gauss.h"
#include "element.h"

struct Element4_2D {
	double **nMatrix, **ksiMatrix, **etaMatrix;
	double dNdX[4][4] = { 0 }, dNdY[4][4] = { 0 };
	double N[4], dKsi[4], dEta[4];
	double *ksi, *eta, *w;
	double detJ = 0;
	int p = 4;

	Element4_2D(int integralPoints = 2);
	~Element4_2D();	
	
	double distance(double x1, double y1, double x2, double y2);
	double* nKsiEta(double ksi, double eta);
	double* dNdEta(double ksi);
	double* dNdKsi(double eta);
	double dXYdEta(int pc, double xy);
	double dXYdEta(double* xy);
	double dXYdKsi(int i, double xy);
	double dXYdKsi(double* xy);
	void fillMatrices();
};