#include "matrix.h"

double* multiply(double* A, double *B, const int N = WIDTH) {
	double* temp = new double[N];
	for (size_t i = 0; i < N; i++)
		temp[i] = A[i] * B[i];

	return temp;
}

double* multiply(double* A, double b, const int N = WIDTH) {
	double* temp = new double[N];
	for (size_t i = 0; i < N; i++)
		temp[i] = A[i] * b;

	return temp;
}

double* add(double* A, double** B, const int M = WIDTH, const int N = HEIGHT) {
	double* temp = new double[M];
	for (size_t i = 0; i < M; i++) {
		for (size_t j = 0; j < N; j++) {
			temp[i] = A[i] + B[i][j];
		}
	}

	return temp;
}

double** add(double** A, double** B, const int M = WIDTH, const int N = HEIGHT) {
	double** temp = new double* [M];
	for (size_t i = 0; i < M; i++) {
		temp[i] = new double[N];
		for (size_t j = 0; j < N; j++) {
			temp[i][j] = A[i][j] + B[i][j];
		}
	}

	return temp;
}

double** divide(double** A, double b, const int M = WIDTH * WIDTH, const int N = HEIGHT * HEIGHT) {
	double** temp = new double* [M];
	for (size_t i = 0; i < M; i++) {
		temp[i] = new double[N];
		for (size_t j = 0; j < N; j++) {
			temp[i][j] = A[i][j] / b;
		}
	}

	return temp;
}

/**
 * Scala macierz r�wna� z macierz� rozwi�za� w jedn�.
 * 
 * @param A - macierz uk�adu r�wna�
 * @param B - macierz rozwi�za�
 * @param n - liczba uk�ad�w r�wna�
 */
double** merge(double** A, double* B, int n) {
	double** AB = new double* [n];

	for (int i = 0; i < n; i++)
		AB[i] = new double[n + 1];
	for (int i = 0; i < n; i++)	
		for (int j = 0; j < n; j++)		
			AB[i][j] = A[i][j];
	
	for (int j = 0; j < n; j++)
		AB[j][n] = B[j];
	
	return AB;
}

/**
 * Zwraca wyznacznik macierzy o wymiarach 2x2.
 *
 * @param matrix - macierz, dla kt�rej wyznaczy� wyznacznik
 * @return (double)wyznacznik macierzy
 */
double detOfMatrix(double matrix[WIDTH][HEIGHT]) {
	return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
}

/**
 * Oblicza odwrotno�� przekazanej macierzy o wymiarach 2x2.
 *
 * @param invM - macierz, do kt�rej zapisa� wynik
 * @param M - macierz do odwr�cenia
 * @param detM - wyznacznik macierzy
 *
 */
void inverseMatrix(double invM[WIDTH][HEIGHT], double M[WIDTH][HEIGHT], double detM) {
	double invDetM = 1 / detM;
	invM[0][0] = M[1][1] * invDetM;
	invM[0][1] = -M[0][1] * invDetM;
	invM[1][0] = -M[1][0] * invDetM;
	invM[1][1] = M[0][0] * invDetM;
}

/**
 * Wype�nia zadan� macierz jej przetransponowan� wersj�.
 *
 * @param src - macierz, z kt�rej ma si� wzorowa�
 * @param dst - macierz, do kt�rej ma zapisa� wynik
 * @param N, M - szeroko�� i wysoko�� macierzy
 */
void transpose(double* src, double* dst, const int M = WIDTH * WIDTH, const int N = HEIGHT * HEIGHT) {
#pragma omp parallel for
	for (int n = 0; n < N * M; n++) {
		int i = n / N;
		int j = n % N;
		dst[n] = src[M * j + i];
	}
}

/*
 * Wy�wietla macierz 16x16.
 *
 * @param M - macierz do wy�wietlenia
 */
void printMatrix(double M[WIDTH * WIDTH * WIDTH * WIDTH][HEIGHT * HEIGHT * HEIGHT * HEIGHT]) {
	for (size_t i = 0; i < 16; i++)
	{
		std::cout << std::setw(8) << M[i][0] << "  " << std::setw(8) << M[i][1] << "  " << std::setw(8) << M[i][2] << "  " << std::setw(8) << M[i][3] << "  ";
		std::cout << std::setw(8) << M[i][4] << "  " << std::setw(8) << M[i][5] << "  " << std::setw(8) << M[i][6] << "  " << std::setw(8) << M[i][7] << "  ";
		std::cout << std::setw(8) << M[i][8] << "  " << std::setw(8) << M[i][9] << "  " << std::setw(8) << M[i][10] << "  " << std::setw(8) << M[i][11] << "  ";
		std::cout << std::setw(8) << M[i][12] << "  " << std::setw(8) << M[i][13] << "  " << std::setw(8) << M[i][14] << "  " << std::setw(8) << M[i][15] << "\n\n";
	}	
}

/*
 * Wy�wietla macierz 4x4.
 *
 * @param M - macierz do wy�wietlenia
 */
void printMatrix(double M[WIDTH * WIDTH][HEIGHT * HEIGHT]) {
	std::cout << std::setw(8) << M[0][0] << "\t" << std::setw(8) << M[0][1] << "\t" << std::setw(8) << M[0][2] << "\t" << std::setw(8) << M[0][3] << "\n";
	std::cout << std::setw(8) << M[1][0] << "\t" << std::setw(8) << M[1][1] << "\t" << std::setw(8) << M[1][2] << "\t" << std::setw(8) << M[1][3] << "\n";
	std::cout << std::setw(8) << M[2][0] << "\t" << std::setw(8) << M[2][1] << "\t" << std::setw(8) << M[2][2] << "\t" << std::setw(8) << M[2][3] << "\n";
	std::cout << std::setw(8) << M[3][0] << "\t" << std::setw(8) << M[3][1] << "\t" << std::setw(8) << M[3][2] << "\t" << std::setw(8) << M[3][3] << "\n\n";
}

/*
 * Wy�wietla macierz 2x2.
 *
 * @param M - macierz do wy�wietlenia
 */
void printMatrix(double M[WIDTH][HEIGHT]) {
	std::cout << std::setw(8) << M[0][0] << "\t" << std::setw(8) << M[0][1] << "\n";
	std::cout << std::setw(8) << M[1][0] << "\t" << std::setw(8) << M[1][1] << "\n\n";
}

/*
 * Wy�wietla macierz 4x1.
 *
 * @param M - macierz do wy�wietlenia
 */
void printMatrix(double M[2 * WIDTH]) {
	std::cout << std::setw(8) << M[0] << "\t" << std::setw(8) << M[1] << "\t" << std::setw(8) << M[2] << "\t" << std::setw(8) << M[3] << "\n\n";
}

/*
 * Wy�wietla macierz MxN.
 *
 * @param matrix - macierz do wy�wietlenia
 * @param M, N - wymiary macierzy
 */
void printMatrix(double **matrix, const int M, const int N) {
	for (size_t i = 0; i < M; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			std::cout << std::setw(10) << matrix[i][j] << "  ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

/*
 * Wy�wietla macierz 1xN.
 *
 * @param matrix - macierz do wy�wietlenia
 * @param N - wymiary macierzy
 */
void printMatrix(double* matrix, const int N) {
	for (size_t i = 0; i < N; i++)	{
		std::cout << std::setw(10) << matrix[i] << "  ";
	}
	std::cout << "\n";
}

/*
 * Wy�wietla 4-kolumnow� macierz z nag��wkami.
 *
 * @param matrix - macierz do wy�wietlenia
 */
void printMatrix4x4(double** matrix) {
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