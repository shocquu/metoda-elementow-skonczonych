#include "matrix.h"

/**
 * Wype³nia zadan¹ macierz jej przetransponowan¹ wersj¹.
 *
 * @param src - macierz, z której ma siê wzorowaæ
 * @param dst - macierz, do której ma zapisaæ wynik
 * @param N, M - szerokoœæ i wysokoœæ macierzy
 */
void transpose(double* src, double* dst, const int N = 4, const int M = 4) {
#pragma omp parallel for
	for (int n = 0; n < N * M; n++) {
		int i = n / N;
		int j = n % N;
		dst[n] = src[M * j + i];
	}
}

/**
 * Zwraca wyznacznik macierzy o wymiarach 2x2.
 *
 * @param matrix - macierz, dla której wyznaczyæ wyznacznik
 * @return (double)wyznacznik macierzy
 */
double detOfMatrix(double matrix[2][2]) {
	return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
}

/**
 * Oblicza odwrotnoœæ przekazanej macierzy o wymiarach 2x2.
 *
 * @param invM - macierz, do której zapisaæ wynik
 * @param M - macierz do odwrócenia
 * @param detM - wyznacznik macierzy
 *
 */
void inverseMatrix(double invM[2][2], double M[2][2], double detM) {
	double invDetM = 1 / detM;
	invM[0][0] = M[1][1] * invDetM;
	invM[0][1] = -M[0][1] * invDetM;
	invM[1][0] = -M[1][0] * invDetM;
	invM[1][1] = M[0][0] * invDetM;
}

/*
 * Wyœwietla macierz 4x4.
 *
 * @param M - macierz do wyœwietlenia
 */
void printMatrix(double M[4][4]) {
	std::cout << std::setw(8) << M[0][0] << "\t" << std::setw(8) << M[0][1] << "\t" << std::setw(8) << M[0][2] << "\t" << std::setw(8) << M[0][3] << "\n";
	std::cout << std::setw(8) << M[1][0] << "\t" << std::setw(8) << M[1][1] << "\t" << std::setw(8) << M[1][2] << "\t" << std::setw(8) << M[1][3] << "\n";
	std::cout << std::setw(8) << M[2][0] << "\t" << std::setw(8) << M[2][1] << "\t" << std::setw(8) << M[2][2] << "\t" << std::setw(8) << M[2][3] << "\n";
	std::cout << std::setw(8) << M[3][0] << "\t" << std::setw(8) << M[3][1] << "\t" << std::setw(8) << M[3][2] << "\t" << std::setw(8) << M[3][3] << "\n\n";
}

/*
 * Wyœwietla macierz 2x2.
 *
 * @param M - macierz do wyœwietlenia
 */
void printMatrix(double M[2][2]) {
	std::cout << std::setw(8) << M[0][0] << "\t" << std::setw(8) << M[0][1] << "\n";
	std::cout << std::setw(8) << M[1][0] << "\t" << std::setw(8) << M[1][1] << "\n\n";
}

/*
 * Wyœwietla macierz 4x1.
 *
 * @param M - macierz do wyœwietlenia
 */
void printMatrix(double M[4]) {
	std::cout << std::setw(8) << M[0] << "\t" << std::setw(8) << M[1] << "\t" << std::setw(8) << M[2] << "\t" << std::setw(8) << M[3] << "\n\n";
}

/*
 * Wyœwietla 4-kolumnow¹ macierz z nag³ówkami.
 *
 * @param matrix - macierz do wyœwietlenia
 */
void printMatrix(double** matrix) {
	std::cout << "+---+-------------+-------------+-------------+-------------+\n";
	std::cout << "| # |           1 |           2 |           3 |           4 |\n";
	std::cout << "+---+-------------+-------------+-------------+-------------+\n";
	for (size_t x = 0; x < 4; x++)
	{
		std::cout << "| " << x + 1 << " ";
		for (size_t y = 0; y < 4; y++)
		{
			std::cout << "|" << std::setw(12) << matrix[x][y] << " ";

		}
		std::cout << "|\n";
		std::cout << "+---+-------------+-------------+-------------+-------------+\n";
	}
	std::cout << "\n";
}