#include "gauss.h"

double Gauss::quadrature1d(double a, double b, vFunctionCall f, int n) {
	double dm = (b - a) / 2;
	double dp = (a + b) / 2;
	double result = 0;

	for (size_t i = 0; i < n + 1; i++) result += A[n - 1][i] * f(dm * x[n - 1][i] + dp);

	return result;
}

double Gauss::quadrature2d(vFunctionCall2 f, int n) {
	double result = 0;

	for (size_t i = 0; i < n + 1; i++) {
		for (size_t j = 0; j < n + 1; j++) result += A[n - 1][i] * A[n - 1][j] * f(x[n - 1][i], x[n - 1][j]);
	}

	return result;
}