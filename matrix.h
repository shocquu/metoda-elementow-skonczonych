#pragma once
#include <iostream>
#include <iomanip>
using namespace std;

#define WIDTH 2
#define HEIGHT 2

double* multiply(double* A, double* B, const int N);
double* multiply(double* A, double b, const int N);
double** add(double** A, double** B, const int M, const int N);
double* add(double* A, double** B, const int M, const int N);
double** divide(double** A, double b, const int M, const int N);
double** merge(double** A, double* B, int n);
double detOfMatrix(double matrix[WIDTH][HEIGHT]);
void inverseMatrix(double invM[WIDTH][HEIGHT], double M[WIDTH][HEIGHT], double detM);
void transpose(double* src, double* dst, const int M, const int N);
void printMatrix(double M[WIDTH * WIDTH * WIDTH * WIDTH][HEIGHT * HEIGHT * HEIGHT * HEIGHT]);
void printMatrix(double M[2 * WIDTH][2 * HEIGHT]);
void printMatrix(double M[WIDTH][HEIGHT]);
void printMatrix(double M[2 * WIDTH]);
void printMatrix(double* matrix, const int N);
void printMatrix(double** matrix, const int M, const int N);
void printMatrix4x4(double** matrix);