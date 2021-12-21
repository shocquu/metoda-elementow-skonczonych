#pragma once
#include <iostream>
#include <iomanip>

#define WIDTH 2
#define WIDTHx4 16
#define HEIGHT 2
#define HEIGHTx4 16

double* multiply(double** A, double* B, const int N = WIDTHx4);
double* multiply(double* A, double* B, const int N = WIDTHx4);
double* multiply(double* A, double b, const int N = WIDTHx4);
double** add(double** A, double** B, const int M = WIDTHx4, const int N = HEIGHTx4);
double* add(double* A, double* B, const int M = WIDTHx4, const int N = HEIGHTx4);
double** divide(double** A, double b, const int M = WIDTHx4, const int N = HEIGHTx4);
double** merge(double** A, double* B, int n = WIDTHx4);
double detOfMatrix(double matrix[WIDTH][HEIGHT]);
void inverseMatrix(double invM[WIDTH][HEIGHT], double M[WIDTH][HEIGHT], double detM);
void transpose(double* src, double* dst, const int M, const int N);
void printMatrix(double M[WIDTHx4][HEIGHTx4]);
void printMatrix(double M[2 * WIDTH][2 * HEIGHT]);
void printMatrix(double M[WIDTH][HEIGHT]);
void printMatrix(double M[2 * WIDTH]);
void printMatrix(double* matrix, const int N = WIDTHx4);
void printMatrix(double** matrix, const int M = WIDTHx4, const int N = HEIGHTx4);
void printMatrix4x4(double** matrix);