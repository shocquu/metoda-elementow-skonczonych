#pragma once
typedef double (*vFunctionCall)(double arg);
typedef double (*vFunctionCall2)(double arg1, double arg2);

class Gauss {
private:
	double x[4][5] = {
		{-0.577350, 0.577350, -1, -1, -1},
		{-0.774597, 0, 0.774597, -1, -1},
		{-0.861136, -0.339981, 0.339981, 0.861136, -1},
		{-0.906180, -0.538469, 0, 0.538469, 0.906180}
	};;
	double A[4][5] = {
		{1, 1, -1, -1, -1},
		{0.555556, 0.888889, 0.555556, -1, -1},
		{0.347855, 0.6521445, 0.6521445, 0.347855, -1},
		{0.236927, 0.478629, 0.568889, 0.478629, 0.236927}
	};;

public:
	double quadrature1d(double a, double b, vFunctionCall f, int n = 2);
	double quadrature2d(vFunctionCall2 f, int n = 2);
};