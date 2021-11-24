#pragma once
typedef double (*vFunctionCall)(double arg);
typedef double (*vFunctionCall2)(double arg1, double arg2);

class Gauss {
private:
	double x[4][5];
	double A[4][5];

public:
	double quadrature1d(double a, double b, vFunctionCall f, int n = 2);
	double quadrature2d(vFunctionCall2 f, int n = 2);
};