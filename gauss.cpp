#include "gauss.h"
#include "matrix.h"

/**
 * Algorytm eliminacji Gauss'a do rozwi¹zywania uk³adów liniowych.
 * 
 * @param AB - scalona macierz uk³adu i rozwi¹zañ
 * @param n - liczba uk³adów równañ
 */
double* Gauss::elimination(double** AB, int n) {
	const double accuracy = 1e-15;
	double* result = new double[n];
	int* vector = new int[n + 1];

	for (int i = 0; i < n + 1; i++)	
		vector[i] = i;
	
	for (int i = 0; i < n - 1; i++) {
		int largest = i;
		bool hasChanged = false;

		for (int j = i + 1; j < n; j++)
			if (fabs(AB[i][vector[largest]]) < fabs(AB[i][vector[j]])) {
				hasChanged = true;
				largest = j;
			}

		if (hasChanged) {
			int pom = vector[i];
			vector[i] = vector[largest];
			vector[largest] = pom;
		}

		for (int j = i + 1; j < n; j++) {
			if (fabs(AB[i][vector[i]]) < accuracy)
				return NULL;
			
			double dzielnik = AB[j][vector[i]] / AB[i][vector[i]];

			for (int k = i + 1; k < n + 1; k++)			
				AB[j][vector[k]] -= (AB[i][vector[k]] * dzielnik);			
		}
	}

	for (int i = n - 1; i >= 0; i--) {
		if (fabs(AB[i][vector[i]]) < accuracy)
			return NULL;
		for (int j = n - 1; j > i; j--)		
			AB[i][n] -= AB[i][vector[j]] * result[vector[j]];
		
		result[vector[i]] = AB[i][n] / AB[i][vector[i]];
	}

	return result;
}

double* Gauss::elimination(double** A, double* B, int n) {
	double** temp = merge(A, B, n);
	return elimination(temp, n);
}

/**
 * Zwraca przybli¿ony wynik ca³kowanej funkcji w przedziale <-1, 1>.
 *
 * @param f - ca³kowana funkcja
 * @param n - liczba wag ?
 */
double Gauss::quadrature1d(vFunctionCall f, int n) {
	double result = 0;

	for (size_t i = 0; i < n + 1; i++)
		result += A[n - 1][i] * f(x[n - 1][i]);

	return result;
}

/**
 * Zwraca przybli¿ony wynik ca³kowanej funkcji w przedziale ze zmian¹ przedzia³u.
 *
 * @param a, b - przedzia³ ca³kowania
 * @param f - ca³kowana funkcja
 * @param n - liczba wag ?
 */
double Gauss::quadrature1d(double a, double b, vFunctionCall f, int n) {
	double dm = (b - a) / 2;
	double dp = (a + b) / 2;
	double result = 0;

	for (size_t i = 0; i < n + 1; i++)
		result += A[n - 1][i] * f(dm * x[n - 1][i] + dp);

	return result;
}

/**
 * Zwraca przybli¿ony wynik ca³kowanej funkcji w przestrzeni 2D.
 *
 * @param f - ca³kowana funkcja
 * @param n - liczba wag ?
 */
double Gauss::quadrature2d(vFunctionCall2 f, int n) {
	double result = 0;

	for (size_t i = 0; i < n + 1; i++)
		for (size_t j = 0; j < n + 1; j++)
			result += A[n - 1][i] * A[n - 1][j] * f(x[n - 1][i], x[n - 1][j]);	

	return result;
}