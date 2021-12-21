#pragma once
#include "grid.h"
#include "element.h"

struct Element4_2D {
	double **nMatrix, **ksiMatrix, **etaMatrix;
	double dNdX[4][4] = { 0 }, dNdY[4][4] = { 0 };
	double N[4], dKsi[4], dEta[4];
	double *ksi, *eta;
	double detJ = 0;
	int p = 4;

	Element4_2D(double* ksi, double* eta, int n);
	~Element4_2D();
	void aggregate(double**& globalH, double**& globalC, double*& globalP, Element currEl);
	void calcHbc(double Hbc[4][4], double P[4], Grid grid, Element currEl, int elIndex, double alpha, double ambientTemp);
	void calcH(double H[4][4], double k);
	void calcC(double C[4][4], double c, double ro, Element currEl, int currElId, Grid grid);
	void calcJ(double J[2][2], int pc, double x, double y);
	void jacobian(Grid grid, Element &currEl);
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